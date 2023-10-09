#include "sys.h"
#include "Primes.h"
#include "utils/ctz.h"
#include "utils/cpu_relax.h"
#include "threadpool/AIThreadPool.h"
#include <tuple>
#include <cmath>
#include <cstring>
#include <thread>

#if defined(__OPTIMIZE__)
#define USE_STOPWATCH 0         // Set to 1 to write execution time information to std::cout.
#endif

#if USE_STOPWATCH
#include "cwds/benchmark.h"
#endif

#ifdef CWDEBUG
#define CHECK_PRIMES 0          // Set to 1 to check that the first 189961811 generated primes are correct.
#define CHECK_SIEVING 0         // Set to 1 to check that bits set in the sieve correspond to integers that are divisible by the current prime.
#else
#define CHECK_PRIMES 0          // Leave at 0.
#define CHECK_SIEVING 0         // idem
#endif

namespace fastprimes {

namespace {

int modular_inverse(int64_t n, int64_t m)
{
  int64_t x0 = 1, x1 = 0;
  int64_t y0 = n, y1 = m;
  while (y1 > 0)
  {
    int64_t q = y0 / y1;
    std::tie(x0, x1) = std::make_tuple(x1, x0 - q * x1);
    std::tie(y0, y1) = std::make_tuple(y1, y0 - q * y1);
  }
  // Ensure the result is positive.
  return (x0 + m) % m;
}

#ifdef CWDEBUG
constexpr int compression_first_prime_second_row = Primes::compression_first_prime + Primes::compression_primorial;
#endif

constexpr std::array<int, compression_repeat> calc_row0()
{
  std::array<int, compression_repeat> row0 = {};
  int col = 0;
  int candidate = Primes::compression_first_prime;
  do
  {
    bool is_divisible = false;
    for (int i = 1; i < Primes::compression; ++i)
      if (candidate % small_primes[i] == 0)
      {
        is_divisible = true;
        break;
      }
    if (!is_divisible)
      row0[col++] = candidate;
    candidate += 2;
  }
  while (col < compression_repeat);
  return row0;
}

constexpr std::array<int, compression_repeat> row0 = calc_row0();

// This is the largest multiplier that is required to keep `offset - row0[col]` positive.
// Since offset = compression_offset_multiplier * prime, that value is maximal when prime
// is as small as possible (prime = compression_first_prime) and row0[col] is as large
// as possible (row0.back()).
constexpr int compression_offset_multiplier = (row0.back() + Primes::compression_first_prime - 1) / Primes::compression_first_prime;

inline prime_t sieve_row_column_to_prime(int row, int column)
{
  return row * static_cast<prime_t>(Primes::compression_primorial) + row0[column];
}

#if CHECK_PRIMES
std::vector<uint32_t> debug_primes;
void debug_init_primes()
{
  std::ifstream ifs("primes_till_4000000000", std::ios::binary);
  if (!ifs.is_open()) {
    DoutFatal(dc::fatal, "Failed to open file primes_till_4000000000 for reading.");
  }

  // Read the size of the vector first.
  size_t size;
  ifs.read(reinterpret_cast<char*>(&size), sizeof(size));

  // Resize the vector and read the data.
  debug_primes.resize(size);
  ifs.read(reinterpret_cast<char*>(debug_primes.data()), size * sizeof(uint32_t));

  ifs.close();
}
#define CHECK_PRIME(prime) \
  ASSERT(debug_pi_ > uint64_t{189961811} || prime == debug_primes[debug_pi_++]);    // We can't test primes larger than 3,999,999,979.
#else
#define CHECK_PRIME(prime) do {} while(0)
#endif // CWDEBUG

#if USE_STOPWATCH
int const cpu = benchmark::Stopwatch::cpu_any;          // The CPU to run Stopwatch on.
double const cpu_frequency = 3612059050.0;              // In cycles per second.
size_t const loopsize = 1000;
size_t const minimum_of = 3;
#endif

std::array<uint16_t, Primes::row0_to_column_size> calc_row0_to_column()
{
  std::array<uint16_t, Primes::row0_to_column_size> result;
#if CW_DEBUG
  std::memset(result.data(), 0, sizeof(result));
#endif
  ASSERT(row0.size() <= std::numeric_limits<uint16_t>::max());
  for (uint16_t i = 0; i < row0.size(); ++i)
  {
    int ri = Primes::n_to_row0_to_column_index(row0[i]);
    ASSERT(result[ri] == 0);
    result[ri] = i;
  }
  return result;
}

} // namespace

// Returns an upper bound on the number of primes that are smaller than or equal to n.
// Just a random function that I thought matched quite well (for n > 1000 or so).
// For n larger than 500,000,000 it is off on average 0.0079% (too large, thus).
//static
integer_t Primes::calc_upper_bound_number_of_primes(integer_t n)
{
  ASSERT(n > 54);       // So that log(n) > 4.
  double logn = std::log(n);
  return std::exp(0.3125 * std::pow(1 / (logn - 4), 1.655) + logn - std::log(logn - 1)) - 4;
}

// COMMENT "precision":
//
// The largest value of `offset - row0[col]` is when `offset` has its largest value
// and row0[col] is at its smallest. The latter happens when col = 0 (at which point
// row0[col] equals compression_first_prime). The former, `offset` at its largest,
// happens when `prime` is at its largest, which is `compression_first_prime_second_row`.
// Note that compression_first_prime_second_row = compression_first_prime + compression_primorial.
//
// Let P = compression_primorial, F = compression_first_prime.
// Then compression_first_prime_second_row = P + F.
// Note that compression_offset_multiplier is more or less (P + F) / F.
//
// And we can write for the largest values involved:
//   prime = P + F,
//   offset = floor((P + F) / F) * (P + F);
//   -row0[col] = -F
//
// The largest possible value of compression_primorial_inverse is prime - 1, or P + F - 1.
//
// Thus the largest possible value of ((offset - row0[col]) * compression_primorial_inverse) is (less than or) equal
//
//   M = (floor((P + F) / F) * (P + F) - F) * (P + F - 1)
//
// which becomes larger than what fits in an int when compression is 6:
//  compression   F     P      M
//  2             5     6      170
//  3             7     30     6408
//  4             11    210    969980
//  5             13    2310   960102882
//  6             17    30030  1595233239472
//
// This means that for compression = 6 we need 64 bit precision when multiplying with the compression_primorial_inverse.
// Calculate the first row that has a multiple of this prime in colum `col`.

template<bool is_row0>
int Primes::wipe_primes_of_row(int row, int number_of_threads)
{
  AIThreadPool& thread_pool = AIThreadPool::instance();
  integer_t const sieve_size = sieve_rows_ * words_per_row;   // The size of the sieve_ in words.
  bool const use_thread_pool = !queue_handle_.undefined();

  sieve_word_t* row_ptr = sieve_ + row;
  int column = 0;
  // Loop over all words in a row, by running over the word index offsets relative to the start of the current row.
  for (unsigned int word_index_offset = 0; word_index_offset < sieve_size; word_index_offset += sieve_rows_)
  {
    // Loop over all bits in a word.
    for (sieve_word_t column_mask = 1; column_mask != 0; column_mask <<= 1, ++column)
    {
      // Did we find the next prime?
      if ((row_ptr[word_index_offset] & column_mask))
      {
        struct WipeWordColData
        {
          prime_t const prime_;
          int const offset_;
          int const compression_primorial_inverse_;
          std::atomic_int col_word_;

          WipeWordColData(prime_t prime) :
            prime_(prime),
            // For row >= 1, prime >= row0[col] so we can simply subtract row0 from prime; this then avoids an overflow.
            // Also see "precision" comment above, for row 0.
            offset_(is_row0 ? compression_offset_multiplier * prime : prime),
            compression_primorial_inverse_(modular_inverse(Primes::compression_primorial, prime)),
            col_word_(0)
          { }
        };

        WipeWordColData data(sieve_row_column_to_prime(row, column));
        ASSERT(is_row0 || data.prime_ > Primes::compression_primorial);

        CHECK_PRIME(data.prime_);

        if constexpr (!is_row0)
        {
          if (data.prime_ > sqrt_max_value_)
          {
            // We're done.
            return column;
          }
        }

        auto wipe_word_column = [this, &data]() -> bool {
          // Run over unprocessed word-columns.
          for (;;)
          {
            // Get a word-column that isn't being wiped yet.
            int col_word = data.col_word_.fetch_add(1, std::memory_order::release);
            // Exit if all word-column have been processed.
            if (col_word >= words_per_row)
              return false;
            unsigned int col_word_offset = col_word * sieve_rows_;
            int col = col_word * sieve_word_bits;
            // Run over all bits in the word(-column).
            for (sieve_word_t col_mask = 1; col_mask != 0; col_mask <<= 1, ++col)
            {
              // Calculate the first row that contains a multiple of data.prime_ in this col(umn).
              uint64_t first_row_with_prime_multiple64 = data.offset_ - row0[col];
              first_row_with_prime_multiple64 *= data.compression_primorial_inverse_;
              int first_row_with_prime_multiple = first_row_with_prime_multiple64 % data.prime_;

              // Run wi (word index) over all words in this word-column that have a bit in this col(umn)
              // that represents an integer that is divisible by prime and reset it.
              for (unsigned int wi = first_row_with_prime_multiple + col_word_offset; wi < sieve_rows_ + col_word_offset; wi += data.prime_)
              {
                sieve_[wi] &= ~col_mask;
#if CHECK_SIEVING
                int debug_row = wi % sieve_rows_;
                int debug_col = (wi / sieve_rows_) * sieve_word_bits + (col % sieve_word_bits);
                prime_t debug_prime = sieve_row_column_to_prime(debug_row, debug_col);
                ASSERT(debug_col == col);
                ASSERT(debug_prime % data.prime_ == 0);
  //            Dout(dc::notice, "Inner loop: setting " << debug_prime << " to 0 because it is " << (debug_prime / data.prime_) << " * " << data.prime_);
#endif
              }
            }
          }
        };
        for (int thread = 0; thread < number_of_threads - 1; ++thread)
        {
          // Get read access to AIThreadPool::m_queues.
          auto queues_access = thread_pool.queues_read_access();
          // Get a reference to one of the queues in m_queues.
          auto& queue = thread_pool.get_queue(queues_access, queue_handle_);
          {
            // Get producer accesses to this queue.
            auto queue_access = queue.producer_access();
            // The queue must be larger than or equal to the number of threads you use.
            ASSERT(queue_access.length() < queue.capacity());
            // Place a lambda in the queue.
            queue_access.move_in(wipe_word_column);
            // Release producer accesses, so another thread can write to this queue again.
          }
          // This function must be called every time move_in was called on a queue that was returned by thread_pool.get_queue.
          queue.notify_one();
          // Release read access to AIThreadPool::m_queues so another thread can use AIThreadPool::new_queue again.
        }
        wipe_word_column();
        // Wait till all threads have finished.
        while (data.col_word_.load(std::memory_order::acquire) < words_per_row + number_of_threads)
          cpu_relax();
        row_ptr[word_index_offset] |= column_mask;
      }
    }
  }
  // Continue with the next row.
  return -1;
}

// Construct a sieve initialized for all primes up till and including max_value.
Primes::Primes(integer_t max_value, int number_of_threads, AIQueueHandle queue_handle) :
  sieve_(nullptr), max_value_(max_value), sqrt_max_value_(std::sqrt(max_value)),
  // Calculate how many integers are not divisible by the first `compression` number of primes that are less than or equal max_value.
  sieve_rows_((max_value - compression_first_prime) / compression_primorial + 1),
  index_(-compression - 1), queue_handle_(queue_handle)
{
  ASSERT(max_value >= compression_first_prime_second_row);

#if USE_STOPWATCH
  benchmark::Stopwatch stopwatch(cpu);          // Declare stopwatch and configure on which CPU it must run.
#endif

#if CHECK_PRIMES
  Debug(debug_init_primes());
  debug_pi_ = compression;
#endif

#if USE_STOPWATCH
  stopwatch.start();
#endif
  // Lets allocate whole rows, even if we only partially need the last one.
  integer_t sieve_size = sieve_rows_ * words_per_row;   // The size of the sieve_ in words.
  index_end_ = sieve_size * sieve_word_bits;            // The size of the sieve_ in bits.

  // Allocate the sieve and fill it with ones.
  sieve_ = (sieve_word_t*)std::malloc(sieve_size * sizeof(sieve_word_t));
  if (!sieve_)
    throw std::bad_alloc();
  // This assumes that unused_bits is zero.
  std::memset(sieve_, 0xff, sieve_size * sizeof(sieve_word_t));
#if USE_STOPWATCH
  stopwatch.stop();
  uint64_t cycles = stopwatch.diff_cycles() - benchmark::Stopwatch::s_stopwatch_overhead;
  float delta = cycles / cpu_frequency;
  std::cout << "Time spent allocating and initializing sieve: " << delta << " seconds." << std::endl;
#endif

  // sieve_ is a one dimensional vector, but can best be throught of as a two dimensional
  // table with width compression_repeat.

  // For example, if compression = 3 then compression_repeat = 8, and sieve_ represents
  // the following integers (all integers not divisible by 2, 3 or 5):
  //
  // row  |col: 0      1      2      3      4      5      6      7
  // -----+--------------------------------------------------------
  //  0 : |  !  7|  @ 11|  # 13|  $ 17|  % 19|  ^ 23|  & 29|  * 31|  <-- row0
  //  1 : |    37|    41|    43|    47|  ! 49|    53|    59|    61|
  //  2 : |    67|    71|    73| @! 77|    79|    83|    89| #! 91|
  //  3 : |    97|   101|   103|   107|   109|   113| $!119|  @121|
  //  4 : |   127|   131| %!133|   137|   139| #@143|   149|   151|
  //  5 : |   157| ^!161|   163|   167|  #169|   173|   179|   181|
  //  6 : | $@187|   191|   193|   197|   199| &!203| %@209|   211|
  //  7 : | *!217| $#221|   223|   227|   229|   233|   239|   241|
  //  8 : | %#247|   251| ^@253|   257|  !259|   263|   269|   271|
  //  9 : |   277|   281|   283|  !287|  $289|   293| ^#299|  !301|
  // 10 : |   307|   311|   313|   317| &@319| %$323|  !329|   331|
  // 11 : |   337| *@341|  !343|   347|   349|   353|   359|  %361|
  // 12 : |   367|  !371|   373| &#377|   379|   383|   389| ^$391|
  // 13 : |   397|   401| *#403|  @407|   409|  !413|   419|   421|
  // 14 : |  !427|   431|   433| ^%437|   439|   443|   449|  @451|
  // 15 : |   457|   461|   463|   467|  !469|  @473|   479|  #481|
  // 16 : |   487|   491| &$493|  !497|   499|   503|   509|  !511|
  // 17 : |  @517|   521|   523| *$527|  ^529|  #533| @!539|   541|
  // 18 : |   547| &%551|  !553|   557|  #559|   563|   569|   571|
  // 19 : |   577|  !581|  @583|   587| *%589|   593|   599|   601|
  // 20 : |   607|  #611|   613|   617|   619|  !623|  $629|   631|
  // 21 : | #!637|   641|   643|   647|  @649|   653|   659|   661|
  // 22 : | &^667|  @671|   673|   677|  !679|   683|  #689|   691|
  // 23 : |  $697|   701|  %703|  !707|   709| *^713|   719|  !721|
  // 24 : |   727|  $731|   733|  @737|   739|   743|  !749|   751|
  // 25 : |   757|   761|  !763|  #767|   769|   773|  %779|  @781|
  // 26 : |   787|  !791|  #793|   797|  $799|  @803|   809|   811|
  // 27 : |  %817|   821|   823|   827|   829| $!833|   839|  &841|
  // 28 : | @!847|  ^851|   853|   857|   859|   863|  @869|  #871|
  // 29 : |   877|   881|   883|   887|  !889|  %893| *&899|  $901|
  // 30 : |   907|   911|  @913|  !917|   919|  #923|   929| %!931|
  // 31 : |   937|   941|  ^943|   947|  #949|   953|  !959|  *961|
  // 32 : | ... etc
  //
  // Here every field contains the integer that is represented by the corresponding bit in sieve_.
  // All such integers divisible by 7 are marked with a '!', those that are divisible by 11 are
  // marked with a '@', and so on. Below we'll say "contains a number" while
  // really it's just a single bit in sieve_ that represents that number.
  //
  // All integers in row 0 have already been calculated during compile time
  // and are stored in row0.
  //
  // In the loop(s) below, we find the next prime (p) and set its corresponding
  // bits to zero. This is done per column; for each column we determine
  // the first row that contains an integer in that column that is divisible
  // by p. For example, for p=13 and col=1, we find that the first number that
  // is divisible by 13 is 221 (= 13*17). While the first number in column 3
  // that is divisible by the prime 11 is 77, row 2.
  //
  // Every subsequent multiple of that prime in the given colum is precisely
  // prime rows lower every time. For example, the vertical distance between
  // every '!' is 7, and between every '$' is 17, etc.
  //
  // Note that the first row does NOT necessarily only contain primes.
  // In this case, the first number that is not a prime is 7^2 = 49 which is
  // larger than 31, the last number in the first row; and thus all numbers in
  // the first row are primes in this case. But that is not the case for
  // compression  = 4 or larger: for compression = 4 the first
  // (top-left) number in the table is 11 and the first non-prime is 11^2 = 121.
  // The first number of row 1 is 11 + (2 * 3 * 5 * 7) = 221. Therefore 121 is
  // part of the first row (and so are 11*13=143, 13*13=169, 11*17=187 and
  // 11*19=209).
  //
  // The first row for a given column that contains a number that is divisible
  // by a prime p is given by: -row0[col] / compression_primorial [MODULO p].
  //
  // For example, we want to know the row in column 2 that contains the first
  // number that is divisible by 19. Note that each time we go one row down,
  // the number is incremented with compression_primorial.
  //
  // row  col=2  MOD 19
  //  0    13      13     = 13 + 0 * 30
  //  1    43       5     = 13 + 1 * 30
  //  2    73      16     = 13 + 2 * 30
  //  3   103       8     = 13 + 3 * 30
  //  4   133       0     = 13 + 4 * 30 --> 4 * 30 = -13 --> 4 = -13 / 30 [MOD 19]
  //
  // The compression_primorial is a constant (only depending on the compression used).
  // But the modulo (p) is not. The row index can be calculated slightly more
  // efficient when p > compression_primorial, which is at least the case for every
  // prime in row 1 or larger (for example, in the case of the above table, with
  // compression = 3, the prime 31 is already larger than the compression_primorial=30,
  // but also 37, the first prime of row 1, is).
  // The first loop deals with these primes (less than compression_primorial) using the
  // fact that row=0, while the second loop makes use of the fact that each prime will
  // be larger than compression_primorial.

  // The sieve_ is a 1D array of words (of type sieve_word_t), having `words_per_row` words per row.
  // For example, for compression = 6 and sieve_word_t = uint64_t, words_per_row = 90 and we
  // can picture the sieve_ as:                        __ compression_repeat = 5760 (number of 'columns').
  //                                                  /
  //                                           55...55
  //                           11  11...11     66...77
  //             00...66  66...22  22...99 ... 99...55
  //     column: 01...23  45...67  89...01     67...89
  //     row=0: [word_00][word_01][word_02]...[word_89]
  //     row=1: [word_90][word_91][word_92]...
  //     row=2:                    010..00 <-- column_mask belonging to column 129.
  //     ...    <--------2------->
  //                      \__ column_word_offset = 2 (belonging to column 129) (called word_index for row 0).
  //

#if USE_STOPWATCH
  stopwatch.start();
#endif
  // Reset all bits in sieve_ corresponding to composite numbers divisible by all primes found in row 0.
  wipe_primes_of_row<true>(0, number_of_threads);
#if USE_STOPWATCH
  stopwatch.stop();
  cycles = stopwatch.diff_cycles() - benchmark::Stopwatch::s_stopwatch_overhead;
  delta = cycles / cpu_frequency;
  std::cout << "Time spent sieving row0: " << delta << " seconds." << std::endl;
#endif

#if USE_STOPWATCH
  stopwatch.start();
#endif
  int row = 0;
  // Loop over all rows, starting with row 1 and find all primes up till
  // and including the first prime that is larger than sqrt_max_value.
  // For all primes less than or equal sqrt_max_value, remove multiples from the sieve_.
  int column;
  while ((column = wipe_primes_of_row<false>(++row, number_of_threads)) == -1)
    ;
#if USE_STOPWATCH
  stopwatch.stop();
  cycles = stopwatch.diff_cycles() - benchmark::Stopwatch::s_stopwatch_overhead;
  delta = cycles / cpu_frequency;
  std::cout << "Time spent sieving row1-" << row << ": " << delta << " seconds." << std::endl;
#endif

#if CHECK_PRIMES
  // Find the word containing the last prime that is less than or equal max_value.
  //
  // For example, with sieve_size = 12 (words)
  //                   sieve_rows_ = 4
  //                   words_per_row = 3
  //
  //    cols: 0-63 | 64-127|128-191
  //         wc:0  | wc:1  | wc:2
  //   r: 0    0       4       8
  //   r: 1    1       5       9
  //   r: 2    2       6      10
  //   r: 3    3       7      11 <-- word_index

  int last_row = sieve_rows_;
  for (;;)
  {
    int last_word_index = sieve_size - sieve_rows_ + --last_row;
    // Run over all columns, right-to-left.
    for (int wc = words_per_row - 1; wc >= 0; --wc)
    {
      int col = wc * sieve_word_bits;
      if (sieve_[last_word_index] != 0 &&
          sieve_row_column_to_prime(last_row, col) <= max_value)
      {
        for (sieve_word_t column_mask = 1; column_mask != 0; column_mask <<= 1, ++col)
        {
          if ((sieve_[last_word_index] & column_mask))
          {
            if (sieve_row_column_to_prime(last_row, col) <= max_value)
              goto found;
            else
              break;
          }
        }
      }
      last_word_index -= sieve_rows_;
    }
  }
found:

#if USE_STOPWATCH
  stopwatch.start();
#endif
  // Add all primes larger than sqrt(max_value) and less than or equal max_value.
  ++column;
  while (row < last_row)
  {
    for (; column < compression_repeat; ++column)
    {
      int column_word_offset = (column / sieve_word_bits) * sieve_rows_;
      sieve_word_t* next_word = sieve_ + row + column_word_offset;
      sieve_word_t column_mask = sieve_word_t{1} << (column % sieve_word_bits);
      if ((*next_word & column_mask))
      {
        prime_t prime = sieve_row_column_to_prime(row, column);
        CHECK_PRIME(prime);
      }
    }
    // Next row.
    ++row;
    column = 0;
  }
#if USE_STOPWATCH
  stopwatch.stop();
  cycles = stopwatch.diff_cycles() - benchmark::Stopwatch::s_stopwatch_overhead;
  delta = cycles / cpu_frequency;
  std::cout << "Time spent adding all primes up till row " << last_row << ": " << delta << " seconds." << std::endl;
#endif

  for (; column < compression_repeat; ++column)
  {
    int column_word_offset = (column / sieve_word_bits) * sieve_rows_;
    sieve_word_t* next_word = sieve_ + row + column_word_offset;
    sieve_word_t column_mask = sieve_word_t{1} << (column % sieve_word_bits);
    if ((*next_word & column_mask))
    {
      prime_t prime = sieve_row_column_to_prime(row, column);
      if (prime > max_value)
        break;
      CHECK_PRIME(prime);
    }
  }
#endif // CHECK_PRIMES
}

prime_t Primes::next_prime()
{
  if (AI_UNLIKELY(++index_ < 0))
    return small_primes[compression + index_];

  // Sieve words, as stored in memory, run top to bottom, left to right.
  // Each "word column" has sieve_rows_ (SR) rows.
  // One row has compression_repeat/sieve_word_bits words (RW) (and compression_repeat columns/bits).
  //
  // column: 0  1  2  3  4  5  6  7    8  9 10 11 12 13 14 15  ...  16 17 18 19 20 21 22 23
  //  row  ---------------------------------------------------------------------------------
  //    0: [         word_0        ] [        word_{SR}      ] ... [ word_{(RW - 1)*SR}    ]
  //    1: [         word_1        ] [        word_{SR+1}    ] ... [ word_{(RW - 1)*SR + 1}]
  //    2: [         word_2        ] [        word_{SR+2}    ] ... [ word_{(RW - 1)*SR + 2}]
  //        ...                                                ...   ...
  //   SR: [         word_{SR-1}   ] [        word_{2*SR-1}  ] ... [ word_{(RW)*SR - 1}    ]
  //
  // index_ runs over all bits of sieve_ from left to right (words left to right, bits inside a word from lsb to msb), top to bottom.
  // For example, if compression_repeat = 24 and sieve_word_bits = 8 we'd have:
  //
  // column: 0  1  2  3  4  5  6  7    8  9 10 11 12 13 14 15  ...  16 17 18 19 20 21 22 23
  //  row  ---------------------------------------------------------------------------------
  //    0: [ 0  1  2  3  4  5  6  7] [ 8  9 10 11 12 13 14 15] ... [16 17 18 19 20 21 22 23]
  //    1: [24 25 26 27 28 29 30 31] [32 33 34 35 36 37 38 39] ... [40 41 42 43 44 45 46 47]
  //    2: [48 49 50 51 52 53 54 55] [56 57 58 59 60 61 62 63] ... [64 65 66 67 68 69 70 71]
  //       ...
  //   SR: [24*{SR-1} ... ] [24*{SR-1}+8 ...] ... [... 24*{SR-1}+23]
  //
  // A single word stores the lowest column (the left-most) in the least significant bit:
  //
  // column: 15 14 13 12 11 10  9  8
  //   bits:  1  1  0  1  0  1  1  1  = 0xd7
  //          0  0  0  1  0  0  0  0  <-- col_mask example (for column = 12).
  //
  // col_mask is word with a single bit set that represents the current column inside the current word.
  // col_word_offset is the word_index corresponding to the current column:
  //
  // column: 0  1  2  3  4  5  6  7    8  9 10 11 12 13 14 15  ...  16 17 18 19 20 21 22 23
  //       [           0           ] [          SR           ] ... [       (RW - 1)*SR     ] <-- col_word_offset
  //
  // The current is therefore given by sieve_[row + col_word_offset].

  static constexpr int words_per_row = compression_repeat / sieve_word_bits;
  int column = index_ % compression_repeat;
  int row    = index_ / compression_repeat;
  sieve_word_t col_mask = sieve_word_t{1} << (column % sieve_word_bits);
  integer_t word_col = column / sieve_word_bits;
  sieve_word_t word = sieve_[row + word_col * sieve_rows_];
  word &= ~(col_mask - 1);      // Remove all bits that represent integers less than the current one.
  // Run over wi (row + word_col * sieve_rows_).
  for (;;)
  {
    if (word)
    {
      column = word_col * sieve_word_bits + utils::ctz(word);
      index_ = int64_t{row} * compression_repeat + column;
      return sieve_row_column_to_prime(row, column);
    }
    if (words_per_row == 1 ||
        AI_UNLIKELY(++word_col == words_per_row))
    {
      if (AI_UNLIKELY(++row == sieve_rows_))
        throw std::out_of_range("Primes: reached end of sieve.");
      if constexpr (words_per_row > 2)
        word_col = 0;
    }
    word = sieve_[row + word_col * sieve_rows_];
  }
}

std::vector<prime_t> Primes::make_vector()
{
  std::vector<prime_t> result;
  result.reserve(calc_upper_bound_number_of_primes(max_value_));
  reset();
  try
  {
    for (;;)
    {
      prime_t p = next_prime();
      if (p > max_value_)
        break;
      result.push_back(p);
    }
  }
  catch (std::out_of_range const&)
  {
  }
  return result;
}

//static
std::array<uint16_t, Primes::row0_to_column_size> const Primes::row0_to_column = calc_row0_to_column();

} // namespace fastprimes
