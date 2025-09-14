#pragma once
// bvn
#if defined(_M_ARM64) || defined(__aarch64__) || defined(_M_ARM) || defined(__arm__)
#include <arm_neon.h>
#elif defined(_M_X64) || defined(__amd64__) || defined(_M_IX86) || defined(__i386__)
#include <immintrin.h>
#else
#error Unsupported platform
#endif
#define bvn 
#include <cstdint>
#include <cstddef>
#include <utility>
#include <cstdlib>
#include <cmath>

#ifdef _KERNEL_MODE
namespace std
{
	template <class _Ty>
	struct remove_reference {
		using type = _Ty;
	};

	template <class _Ty>
	struct remove_reference<_Ty&> {
		using type = _Ty;
	};

	template <class _Ty>
	struct remove_reference<_Ty&&> {
		using type = _Ty;
	};

	template <class _Ty>
	using remove_reference_t = typename remove_reference<_Ty>::type;

	template <class _Ty>
	struct remove_const {
		using type = _Ty;
	};

	template <class _Ty>
	struct remove_const<const _Ty> {
		using type = _Ty;
	};

	template <class _Ty>
	using remove_const_t = typename remove_const<_Ty>::type;

	template<class T>
	struct make_unsigned;

	template<>
	struct make_unsigned<char> { using type = unsigned char; };
	
	template<>
	struct make_unsigned<signed char> { using type = unsigned char; };
	
	template<>
	struct make_unsigned<unsigned char> { using type = unsigned char; };

	template<class T, T... Ints>
	struct integer_sequence {
		using value_type = T;
		static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
	};

	template<std::size_t... Ints>
	using index_sequence = integer_sequence<std::size_t, Ints...>;

	template<class T, std::size_t N>
	struct make_integer_sequence_impl;

	template<std::size_t N>
	using make_index_sequence = typename make_integer_sequence_impl<std::size_t, N>::type;

	template<class T, T N>
	using make_integer_sequence = typename make_integer_sequence_impl<T, N>::type;
}
#else
#include <type_traits>
#endif

#ifdef _MSC_VER
#define ENCRYP_FORCEINLINE __forceinline
#else
#define ENCRYP_FORCEINLINE __attribute__((always_inline)) inline
#endif

#define encryp(x) do { \
    static volatile int seed = 10003; \
    \
     \
    if (rand() % 2) { \
         \
        for (int _junk_iter = 0; _junk_iter < 6; _junk_iter++) { \
            volatile int var1 = 0xA024093; \
            volatile float var2 = 3.14159265358979323846f; \
            volatile double var3 = 2.718281828459045; \
            volatile long long var4 = 0xA59950938903403LL; \
            volatile short arr[128]; \
            volatile char char_arr[256]; \
            volatile int jmp_table[32]; \
            volatile int jmp_history[64] = {0}; \
            volatile int jmp_index = 0; \
            volatile int layer = 0; \
            volatile int dummy1 = 0x589403; \
            volatile float dummy2 = 1.234567f; \
            volatile double dummy3 = 9.876543210; \
            volatile long long dummy4 = 0x4835804; \
            volatile short dummy_arr1[64]; \
            volatile char dummy_arr2[128]; \
            volatile int dummy_arr3[32]; \
            \
            for (int i = 0; i < 32; i++) { \
                jmp_table[i] = (seed ^ i) % 32; \
                dummy_arr3[i] = (dummy1 ^ i) * 0x9E3779B9; \
            } \
            seed = (seed * 0x8088405 + 1) & 0xFFFFFFFF; \
            dummy1 = (dummy1 * 0x41C64E6D + 0x3039) & 0xFFFFFFFF; \
            \
            for (volatile int _a_ = 0; _a_ < 100; ++_a_) { \
                var1 ^= (_a_ * 0x1337); \
                dummy1 ^= (_a_ * 0x5EED); \
                if (_a_ % 11 == 0 && layer < 3) { \
                    jmp_history[jmp_index++ % 64] = _a_; \
                    layer++; \
                } \
                \
                for (volatile int _b_ = 0; _b_ < 30; ++_b_) { \
                    var2 *= (1.0f + (_b_ * 0.01f)); \
                    dummy2 *= (1.0f + (_b_ * 0.005f)); \
                    if (_b_ % 15 == 1 && layer > 0) { \
                        jmp_history[jmp_index++ % 64] = _b_; \
                        layer--; \
                    } \
                    \
                    if (((var1 ^ _a_) & (_b_ + 1)) % 7 == 0) { \
                        var3 += var2 / (1.0 + _a_); \
                        dummy3 += dummy2 / (1.0 + _a_); \
                        if (((var1 + _a_ * _b_) % (_b_ + 5)) == 0) { \
                            var4 ^= (0x49058309LL << (_a_ % 16)); \
                            dummy4 ^= (0x13579BDFLL << (_b_ % 16)); \
                            switch ((var1 ^ (_a_ * _b_)) % 20) { \
                                case 0: arr[_a_ % 128] = _b_; dummy_arr1[_a_ % 64] = _b_; break; \
                                case 1: var1 = ~var1; dummy1 = ~dummy1; break; \
                                case 2: var2 = -var2; dummy2 = -dummy2; break; \
                                case 3: var3 *= 0.5; dummy3 *= 0.7; break; \
                                case 4: var4 >>= 1; dummy4 <<= 1; break; \
                                case 5: char_arr[(_a_ + _b_) % 256] = _a_ ^ _b_; dummy_arr2[(_a_ + _b_) % 128] = _a_ + _b_; break; \
                                case 6: var1 = var1 | (1 << (_a_ % 32)); dummy1 = dummy1 & ~(1 << (_b_ % 32)); break; \
                                case 7: var2 += var3 / 1000.0f; dummy2 -= dummy3 / 2000.0f; break; \
                                case 8: var3 = (var3 > 1000) ? 0 : var3 * 2; dummy3 = (dummy3 < -500) ? 1000 : dummy3 / 2; break; \
                                case 9: var4 = (var4 * 7) % 0xFFFFFFFFFFFFLL; dummy4 = (dummy4 * 11) % 0xFFFFFFFFFFFFLL; break; \
                                case 10: arr[(_a_ * _b_) % 128] = _a_ + _b_; dummy_arr1[(_a_ + _b_) % 64] = _a_ * _b_; break; \
                                case 11: var1 = var1 & ~(1 << (_b_ % 32)); dummy1 = dummy1 | (1 << (_a_ % 32)); break; \
                                case 12: var2 *= (_a_ % 2) ? 1.5f : 0.5f; dummy2 *= (_b_ % 2) ? 0.8f : 1.2f; break; \
                                case 13: var3 += sin((double)_a_ / (double)(_b_ + 1)); dummy3 += cos((double)_b_ / (double)(_a_ + 1)); break; \
                                case 14: var4 ^= (0x45384LL << (_b_ % 8)); dummy4 ^= (0x97531LL << (_a_ % 8)); break; \
                                case 15: arr[(_a_ + _b_ * 3) % 128] = _a_ * _b_; dummy_arr1[(_b_ + _a_ * 2) % 64] = _a_ - _b_; break; \
                                case 16: jmp_history[jmp_index++ % 64] = 16; dummy_arr3[(_a_ + _b_) % 32] = _a_ ^ _b_; break; \
                                case 17: var1 = var1 ^ (seed * _a_ * _b_); dummy1 = dummy1 ^ (seed * _b_ * _a_); break; \
                                case 18: if (layer < 3) { layer++; } dummy4 += _a_ * _b_; break; \
                                case 19: arr[(_a_ * _b_) % 128] ^= 0xFFFF; dummy_arr2[(_a_ + _b_) % 128] ^= 0xAA; break; \
                            } \
                        } \
                    } \
                    \
                    if ((_a_ ^ _b_) % 3 == 0) { var1 += _a_ * _b_; dummy1 -= _a_ + _b_; } \
                    if ((_a_ + _b_) % 5 == 0) { var3 *= 1.001; dummy3 *= 0.999; } \
                    if ((_a_ * _b_) % 7 == 0) { var4 ^= (1ULL << (_a_ % 63)); dummy4 ^= (1ULL << (_b_ % 63)); } \
                    \
                    for (volatile int _c_ = 0; _c_ < 5 && _c_ < _b_; ++_c_) { \
                        var1 = (var1 * 0x17489 + 0x24A63) & 0xFFFFFFFF; \
                        dummy1 = (dummy1 * 0x25691 + 0x13B57) & 0xFFFFFFFF; \
                        arr[(_a_ + _b_ + _c_) % 128] = _c_ ^ _a_ ^ _b_; \
                        dummy_arr1[(_a_ + _c_) % 64] = _c_ + _a_ + _b_; \
                        if (_c_ % 3 == 0) { \
                            var2 *= exp(sin((float)_c_ / 10.0f)); \
                            dummy2 *= log(cos((float)_c_ / 5.0f) + 1.0f); \
                        } else if (_c_ % 3 == 1) { \
                            var3 = tan(var3 / 100.0) * 10.0; \
                            dummy3 = sinh(dummy3 / 200.0) * 20.0; \
                        } else { \
                            var4 ^= (seed ^ 0x0980348509LL) * _c_; \
                            dummy4 ^= (seed ^ 0x1234567890LL) * (_c_ + 1); \
                        } \
                    } \
                } \
            } \
            \
            if (var1 == 0x12345678 && var4 == 0x87654321LL) { \
                for (int i = 0; i < 128; i++) { arr[i] = 0; dummy_arr1[i % 64] = 0; } \
            } \
        } \
    } else { \
         \
        for (int _alt_iter = 0; _alt_iter < 20; _alt_iter++) { \
            volatile int var1 = 0x3985; \
            volatile float var2 = 3.14159265358979323846f; \
            volatile double var3 = 2.718281828459045; \
            volatile long long var4 = 0x49534834804LL; \
            volatile short arr[128]; \
            volatile char char_arr[256]; \
            volatile int jmp_table[32]; \
            volatile int jmp_history[64] = {0}; \
            volatile int jmp_index = 0; \
            volatile int layer = 0; \
            volatile int extra1 = 0x48392; \
            volatile float extra2 = 2.71828f; \
            volatile double extra3 = 1.41421356237; \
            volatile long long extra4 = 0x1122334455667788LL; \
            volatile char extra_arr1[96]; \
            volatile int extra_arr2[48]; \
            volatile short extra_arr3[80]; \
            \
            for (int i = 0; i < 32; i++) { \
                jmp_table[i] = (seed ^ i) % 32; \
                extra_arr2[i % 48] = (extra1 ^ (i * 3)) + seed; \
            } \
            seed = (seed * 0x8088405 + 1) & 0xFFFFFFFF; \
            extra1 = (extra1 * 0x6AC690C5 + 0x5F3759DF) & 0xFFFFFFFF; \
            \
            for (volatile int _a_ = 0; _a_ < 50; ++_a_) { \
                var1 ^= (_a_ * 0x1337); \
                extra1 ^= (_a_ * 0x2468); \
                if (_a_ % 11 == 0 && layer < 3) { \
                    jmp_history[jmp_index++ % 64] = _a_; \
                    layer++; \
                } \
                \
                for (volatile int _b_ = 0; _b_ < 15; ++_b_) { \
                    var2 *= (1.0f + (_b_ * 0.01f)); \
                    extra2 *= (1.0f + (_b_ * 0.007f)); \
                    if (_b_ % 15 == 1 && layer > 0) { \
                        jmp_history[jmp_index++ % 64] = _b_; \
                        layer--; \
                    } \
                    \
                    if (((var1 ^ _a_) & (_b_ + 1)) % 7 == 0) { \
                        var3 += var2 / (1.0 + _a_); \
                        extra3 += extra2 / (2.0 + _b_); \
                        if (((var1 + _a_ * _b_) % (_b_ + 5)) == 0) { \
                            var4 ^= (0x85498340430LL << (_a_ % 16)); \
                            extra4 ^= (0x24681357ACEULL << (_b_ % 16)); \
                            switch ((var1 ^ (_a_ * _b_)) % 10) { \
                                case 0: arr[_a_ % 128] = _b_; extra_arr3[_a_ % 80] = _b_ * 2; break; \
                                case 1: var1 = ~var1; extra1 = ~extra1; break; \
                                case 2: var2 = -var2; extra2 = -extra2; break; \
                                case 3: var3 *= 0.5; extra3 *= 0.6; break; \
                                case 4: var4 >>= 1; extra4 <<= 2; break; \
                                case 5: char_arr[(_a_ + _b_) % 256] = _a_ ^ _b_; extra_arr1[(_a_ + _b_) % 96] = _a_ + _b_; break; \
                                case 6: var1 = var1 | (1 << (_a_ % 32)); extra1 = extra1 & ~(1 << (_b_ % 32)); break; \
                                case 7: var2 += var3 / 1000.0f; extra2 -= extra3 / 1500.0f; break; \
                                case 8: var3 = (var3 > 1000) ? 0 : var3 * 2; extra3 = (extra3 < -300) ? 500 : extra3 / 3; break; \
                                case 9: var4 = (var4 * 7) % 0xFFFFFFFFFFFFLL; extra4 = (extra4 * 13) % 0xFFFFFFFFFFFFLL; break; \
                            } \
                        } \
                    } \
                } \
            } \
        } \
    } \
    \
     \
    for (int _bug_iter = 0; _bug_iter < 14; _bug_iter++) { \
        if (_bug_iter % 2 == 0) { \
             \
            volatile int var1 = 0xA4385834; \
            volatile float var2 = 3.14159265358979323846f; \
            volatile double var3 = 2.718281828459045; \
            volatile long long var4 = 0xA5893485039503LL; \
            volatile short arr[64]; \
            volatile char char_arr[128]; \
            volatile int misc1 = 0x87654321; \
            volatile float misc2 = 1.618033988f; \
            volatile double misc3 = 0.577215664901; \
            volatile long long misc4 = 0xFEDCBA9876543210LL; \
            volatile char misc_arr1[112]; \
            volatile short misc_arr2[56]; \
            volatile int misc_arr3[28]; \
            \
            for (volatile int _a_ = 0; _a_ < 50; ++_a_) { \
                var1 ^= (_a_ * 0x1337); \
                misc1 ^= (_a_ * 0x9ABC); \
                for (volatile int _b_ = 0; _b_ < 15; ++_b_) { \
                    var2 *= (1.0f + (_b_ * 0.01f)); \
                    misc2 *= (1.0f + (_b_ * 0.013f)); \
                    if (((var1 ^ _a_) & (_b_ + 1)) % 7 == 0) { \
                        var3 += var2 / (1.0 + _a_); \
                        misc3 += misc2 / (3.0 + _a_); \
                        if (((var1 + _a_ * _b_) % (_b_ + 5)) == 0) { \
                            var4 ^= (0x5945L << (_a_ % 16)); \
                            misc4 ^= (0x9383L << (_b_ % 16)); \
                            arr[_a_ % 64] = _b_; \
                            char_arr[(_a_ + _b_) % 128] = _a_ ^ _b_; \
                            misc_arr1[(_a_ + _b_) % 112] = _a_ - _b_; \
                            misc_arr2[(_a_ * _b_) % 56] = _a_ + _b_; \
                            misc_arr3[(_a_ ^ _b_) % 28] = _a_ * _b_; \
                        } \
                    } \
                } \
            } \
        } else { \
             \
            volatile int junk_var1 = 0x93485030; \
            volatile double junk_var2 = 1.414213562373095; \
            volatile long long junk_var3 = 0x90453803495043LL; \
            volatile float junk_var4 = 0.693147180f; \
            volatile int junk_var5 = 0x55AA55AA; \
            volatile short junk_var6 = 0x1234; \
            volatile char junk_arr1[144]; \
            volatile int junk_arr2[72]; \
            volatile short junk_arr3[88]; \
            volatile long long junk_arr4[36]; \
            \
            for (volatile int _x_ = 0; _x_ < 30; ++_x_) { \
                junk_var1 ^= (_x_ * 0x4893); \
                junk_var5 ^= (_x_ * 0x6789); \
                for (volatile int _y_ = 0; _y_ < 10; ++_y_) { \
                    junk_var2 *= (1.0 + (_y_ * 0.02)); \
                    junk_var4 *= (1.0f + (_y_ * 0.015f)); \
                    if ((junk_var1 + _x_ * _y_) % 7 == 0) { \
                        junk_var3 ^= (0x3495830LL << (_x_ % 8)); \
                        junk_arr4[(_x_ + _y_) % 36] = junk_var3 + _x_ * _y_; \
                    } \
                    if ((_x_ + _y_) % 3 == 0) { \
                        junk_arr1[(_x_ * _y_) % 144] = _x_ ^ _y_; \
                        junk_arr2[(_x_ + _y_) % 72] = junk_var1 + junk_var5; \
                        junk_arr3[(_x_ - _y_ + 50) % 88] = junk_var6 * _x_; \
                    } \
                } \
            } \
        } \
    } \
} while(0)

namespace encryp_core
{
	template<class _Ty>
	using clean_type = typename std::remove_const_t<std::remove_reference_t<_Ty>>;

	namespace detail {
		template<std::size_t Size>
		ENCRYP_FORCEINLINE constexpr std::size_t buffer_size()
		{
			return ((Size / 16) + (Size % 16 != 0)) * 2;
		}

		template<std::uint32_t Seed>
		ENCRYP_FORCEINLINE constexpr std::uint32_t encryp_key4() noexcept
		{
			std::uint32_t value = Seed;
			for(char c : __TIME__)
				value = static_cast<std::uint32_t>((value ^ c) * 16777619ull);
			return value;
		}

		template<std::size_t S>
		ENCRYP_FORCEINLINE constexpr std::uint64_t encryp_key8()
		{
			constexpr auto first_part = encryp_key4<2166136261 + S>();
			constexpr auto second_part = encryp_key4<first_part>();
			return (static_cast<std::uint64_t>(first_part) << 32) | second_part;
		}

		template<std::size_t N, class CharT>
		ENCRYP_FORCEINLINE constexpr std::uint64_t
		load_encryped_str8(std::uint64_t key, std::size_t idx, const CharT* str) noexcept
		{
			using cast_type = typename std::make_unsigned<CharT>::type;
			constexpr auto value_size = sizeof(CharT);
			constexpr auto idx_offset = 8 / value_size;

			std::uint64_t value = key;
			for(std::size_t i = 0; i < idx_offset && i + idx * idx_offset < N; ++i)
				value ^=
					(std::uint64_t{ static_cast<cast_type>(str[i + idx * idx_offset]) }
					 << ((i % idx_offset) * 8 * value_size));

			return value;
		}

		ENCRYP_FORCEINLINE std::uint64_t load_from_reg(std::uint64_t value) noexcept
		{
#if defined(__clang__) || defined(__GNUC__)
			asm("" : "=r"(value) : "0"(value) :);
			return value;
#else
			volatile std::uint64_t reg = value;
			return reg;
#endif
		}
	}

	template <int _size, char _key1, char _key2, typename T>
	class encryp_basic
	{
	public:
		__forceinline constexpr encryp_basic(T* data)
		{
			encryp_data(data);
		}

		__forceinline T* get()
		{
			return _storage;
		}

		__forceinline int size()
		{
			return _size;
		}

		__forceinline char key()
		{
			return _key1;
		}

		__forceinline T* (encryp)()
		{
			if (!isEncryped())
				encryp_data(_storage);
			return _storage;
		}

		__forceinline T* decryp()
		{
			if (isEncryped())
				encryp_data(_storage);
			return _storage;
		}

		__forceinline bool isEncryped()
		{
			return _storage[_size - 1] != 0;
		}

		__forceinline void clear()
		{
			for (int i = 0; i < _size; i++)
			{
				_storage[i] = 0;
			}
		}

		__forceinline operator T* ()
		{
			decryp();
			return _storage;
		}

	private:
		__forceinline constexpr void encryp_data(T* data)
		{
			for (int i = 0; i < _size; i++)
			{
				_storage[i] = data[i] ^ (_key1 + i % (1 + _key2));
			}
		}

		T _storage[_size]{};
	};

	template<class CharT, std::size_t Size, class Keys, class Indices>
	class encryp_advanced;

	template<class CharT, std::size_t Size, std::uint64_t... Keys, std::size_t... Indices>
	class encryp_advanced<CharT, Size, std::integer_sequence<std::uint64_t, Keys...>, std::index_sequence<Indices...>> {
#ifndef ENCRYP_DISABLE_AVX_INTRINSICS
		constexpr static inline std::uint64_t alignment = ((Size > 16) ? 32 : 16);
#else
		constexpr static inline std::uint64_t alignment = 16;
#endif

		alignas(alignment) std::uint64_t _storage[sizeof...(Keys)];

	public:
		using value_type = CharT;
		using size_type = std::size_t;
		using pointer = CharT*;
		using const_pointer = const CharT*;

		template<class L>
		ENCRYP_FORCEINLINE encryp_advanced(L l, std::integral_constant<std::size_t, Size>, std::index_sequence<Indices...>) noexcept
			: _storage{ detail::load_from_reg((std::integral_constant<std::uint64_t, detail::load_encryped_str8<Size>(Keys, Indices, l())>::value))... }
		{}

		ENCRYP_FORCEINLINE constexpr size_type size() const noexcept
		{
			return Size - 1;
		}

		ENCRYP_FORCEINLINE void (encryp)() noexcept
		{
#if defined(__clang__)
			alignas(alignment)
				std::uint64_t arr[]{ detail::load_from_reg(Keys)... };
			std::uint64_t* keys =
				(std::uint64_t*)detail::load_from_reg((std::uint64_t)arr);
#else
			alignas(alignment) std::uint64_t keys[]{ detail::load_from_reg(Keys)... };
#endif

#if defined(_M_ARM64) || defined(__aarch64__) || defined(_M_ARM) || defined(__arm__)
#if defined(__clang__)
			((Indices >= sizeof(_storage) / 16 ? static_cast<void>(0) : __builtin_neon_vst1q_v(
				reinterpret_cast<uint64_t*>(_storage) + Indices * 2,
				veorq_u64(__builtin_neon_vld1q_v(reinterpret_cast<const uint64_t*>(_storage) + Indices * 2, 51),
					__builtin_neon_vld1q_v(reinterpret_cast<const uint64_t*>(keys) + Indices * 2, 51)),
				51)), ...);
#else
			((Indices >= sizeof(_storage) / 16 ? static_cast<void>(0) : vst1q_u64(
				reinterpret_cast<uint64_t*>(_storage) + Indices * 2,
				veorq_u64(vld1q_u64(reinterpret_cast<const uint64_t*>(_storage) + Indices * 2),
					vld1q_u64(reinterpret_cast<const uint64_t*>(keys) + Indices * 2)))), ...);
#endif
#elif !defined(ENCRYP_DISABLE_AVX_INTRINSICS)
			((Indices >= sizeof(_storage) / 32 ? static_cast<void>(0) : _mm256_store_si256(
				reinterpret_cast<__m256i*>(_storage) + Indices,
				_mm256_xor_si256(
					_mm256_load_si256(reinterpret_cast<const __m256i*>(_storage) + Indices),
					_mm256_load_si256(reinterpret_cast<const __m256i*>(keys) + Indices)))), ...);

			if constexpr(sizeof(_storage) % 32 != 0)
				_mm_store_si128(
					reinterpret_cast<__m128i*>(_storage + sizeof...(Keys) - 2),
					_mm_xor_si128(_mm_load_si128(reinterpret_cast<const __m128i*>(_storage + sizeof...(Keys) - 2)),
						_mm_load_si128(reinterpret_cast<const __m128i*>(keys + sizeof...(Keys) - 2))));
#else
			((Indices >= sizeof(_storage) / 16 ? static_cast<void>(0) : _mm_store_si128(
				reinterpret_cast<__m128i*>(_storage) + Indices,
				_mm_xor_si128(_mm_load_si128(reinterpret_cast<const __m128i*>(_storage) + Indices),
					_mm_load_si128(reinterpret_cast<const __m128i*>(keys) + Indices)))), ...);
#endif
		}

		ENCRYP_FORCEINLINE const_pointer get() const noexcept
		{
			return reinterpret_cast<const_pointer>(_storage);
		}

		ENCRYP_FORCEINLINE pointer get() noexcept
		{
			return reinterpret_cast<pointer>(_storage);
		}

		ENCRYP_FORCEINLINE pointer get_encryped() noexcept
		{
#if defined(__clang__)
			alignas(alignment)
				std::uint64_t arr[]{ detail::load_from_reg(Keys)... };
			std::uint64_t* keys =
				(std::uint64_t*)detail::load_from_reg((std::uint64_t)arr);
#else
			alignas(alignment) std::uint64_t keys[]{ detail::load_from_reg(Keys)... };
#endif

#if defined(_M_ARM64) || defined(__aarch64__) || defined(_M_ARM) || defined(__arm__)
#if defined(__clang__)
			((Indices >= sizeof(_storage) / 16 ? static_cast<void>(0) : __builtin_neon_vst1q_v(
				reinterpret_cast<uint64_t*>(_storage) + Indices * 2,
				veorq_u64(__builtin_neon_vld1q_v(reinterpret_cast<const uint64_t*>(_storage) + Indices * 2, 51),
					__builtin_neon_vld1q_v(reinterpret_cast<const uint64_t*>(keys) + Indices * 2, 51)),
				51)), ...);
#else
			((Indices >= sizeof(_storage) / 16 ? static_cast<void>(0) : vst1q_u64(
				reinterpret_cast<uint64_t*>(_storage) + Indices * 2,
				veorq_u64(vld1q_u64(reinterpret_cast<const uint64_t*>(_storage) + Indices * 2),
					vld1q_u64(reinterpret_cast<const uint64_t*>(keys) + Indices * 2)))), ...);
#endif
#elif !defined(ENCRYP_DISABLE_AVX_INTRINSICS)
			((Indices >= sizeof(_storage) / 32 ? static_cast<void>(0) : _mm256_store_si256(
				reinterpret_cast<__m256i*>(_storage) + Indices,
				_mm256_xor_si256(
					_mm256_load_si256(reinterpret_cast<const __m256i*>(_storage) + Indices),
					_mm256_load_si256(reinterpret_cast<const __m256i*>(keys) + Indices)))), ...);

			if constexpr(sizeof(_storage) % 32 != 0)
				_mm_store_si128(
					reinterpret_cast<__m128i*>(_storage + sizeof...(Keys) - 2),
					_mm_xor_si128(_mm_load_si128(reinterpret_cast<const __m128i*>(_storage + sizeof...(Keys) - 2)),
						_mm_load_si128(reinterpret_cast<const __m128i*>(keys + sizeof...(Keys) - 2))));
#else
			((Indices >= sizeof(_storage) / 16 ? static_cast<void>(0) : _mm_store_si128(
				reinterpret_cast<__m128i*>(_storage) + Indices,
				_mm_xor_si128(_mm_load_si128(reinterpret_cast<const __m128i*>(_storage) + Indices),
					_mm_load_si128(reinterpret_cast<const __m128i*>(keys) + Indices)))), ...);
#endif

			return (pointer)(_storage);
		}
	};

	template<class L, std::size_t Size, std::size_t... Indices>
	encryp_advanced(L l, std::integral_constant<std::size_t, Size>, std::index_sequence<Indices...>) -> encryp_advanced<
		std::remove_const_t<std::remove_reference_t<decltype(l()[0])>>,
		Size,
		std::integer_sequence<std::uint64_t, detail::encryp_key8<Indices>()...>,
		std::index_sequence<Indices...>>;

	template<class StrType>
	ENCRYP_FORCEINLINE auto make_encryp_advanced(StrType&& str)
	{
		return encryp_advanced([]() { return str; }, 
			std::integral_constant<std::size_t, sizeof(str) / sizeof(*str)>{}, 
			std::make_index_sequence<detail::buffer_size<sizeof(str)>()>{});
	}

	template<class StrType>
	ENCRYP_FORCEINLINE auto make_encryp_basic(StrType&& str)
	{
		return encryp_basic<sizeof(str) / sizeof(str[0]), __TIME__[4], __TIME__[7], 
			clean_type<decltype(str[0])>>((clean_type<decltype(str[0])>*)str);
	}
}

#define encrypt(str) []() { \
	static auto encryped = encryp_core::make_encryp_basic(str); \
	return encryped; }()