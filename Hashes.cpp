#include "Hashes.h"

#include "Random.h"

#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
#include <assert.h>
//#include <emmintrin.h>
//#include <xmmintrin.h>


#include <stdint.h>
#include <stdio.h>

// Rotate left (cyclic shift) within 32 bits.
static inline uint32_t rotateLeft32(uint32_t x, unsigned r) {
  return (x << r) | (x >> (32 - r));
}

// Taking the four bytes, computing their partial sums (casted to mod 256 (8-bit)), then
// XOR'ing them with 8-bit shifts to "append" them to a 32-bit int.
static inline uint32_t packPartialSums(const uint8_t b0,
                                     const uint8_t b1,
                                     const uint8_t b2,
                                     const uint8_t b3) {
  // c0 = b0
  // c1 = b0 + b1 (mod 256)
  // c2 = b0 + b1 + b2 (mod 256)
  // c3 = b0 + b1 + b2 + b3 (mod 256)
  uint8_t c0 = b0;
  uint8_t c1 = static_cast<uint8_t>(b0 + b1);
  uint8_t c2 = static_cast<uint8_t>(b0 + b1 + b2);
  uint8_t c3 = static_cast<uint8_t>(b0 + b1 + b2 + b3);

  // Taking a blank 32-bit integer.
  uint32_t out = 0;

  // c0 goes into the least significant byte.
  out ^= (uint32_t)c0;

  // c1 goes into the second byte from the shift, c2 goes into the third byte, etc.
  out ^= (uint32_t)c1 << 8;
  out ^= (uint32_t)c2 << 16;
  out ^= (uint32_t)c3 << 24;
  return out;
}

// Our 32-bit toy hash with cumulative partial sums and a cyclic shift.
uint32_t StudentHash(std::istream &in) {
  // Initialization vector.
  uint32_t hashState = 0xdeadbeef;
  
  // This is going to be looped endlessly until we reach the end of the message.
  while (true) {
      // Read up to four bytes of the buffer, and we will only do four at a time.
      char buffer[4] = {0, 0, 0, 0};
      in.read(buffer, 4);

      // This is a counter to track how much buffer is left to be proccessed.
      std::streamsize bytesRead = in.gcount();
      if (bytesRead == 0) {
          // No more data.
          break;
      }
      
      // For the four values in the buffer, I will be padding with zeros in case the
      // buffer is less than four. So that the partial sums can still be computed.
      uint8_t b0 = (bytesRead > 0) ? (uint8_t)buffer[0] : 0;
      uint8_t b1 = (bytesRead > 1) ? (uint8_t)buffer[1] : 0;
      uint8_t b2 = (bytesRead > 2) ? (uint8_t)buffer[2] : 0;
      uint8_t b3 = (bytesRead > 3) ? (uint8_t)buffer[3] : 0;
      
      // Compute cumulative partial sums mod 256 and pack them into a 32-bit block.
      uint32_t sumBlock = packPartialSums(b0, b1, b2, b3);
      
      // --- Combine & Mix Steps ---
      
      // 1) Rotate hashState left by a constant (e.g., 5 bits).
      hashState = rotateLeft32(hashState, 5);
      
      // 2) XOR the cumulative-sum block into the state.
      hashState ^= sumBlock;
      
      // 3) Multiply by a prime constant mod 2^32 to further mix the bits.
      hashState *= 0x1000193;
      
      // If we did not read a full 4-byte block, we're finished.
      if (bytesRead < 4) {
          break;
      }
  }
  
    // Write final 32-bit digest
    *(uint32_t*)out = hashState;
}

//----------------------------------------------------------------------------
// fake / bad hashes

void BadHash ( const void * key, int len, uint32_t seed, void * out )
{
  uint32_t h = seed;

  const uint8_t * data = (const uint8_t*)key;

  for(int i = 0; i < len; i++)
  {
    h ^= h >> 3;
    h ^= h << 5;
    h ^= data[i];
  }

  *(uint32_t*)out = h;
}

void sumhash ( const void * key, int len, uint32_t seed, void * out )
{
  uint32_t h = seed;

  const uint8_t * data = (const uint8_t*)key;

  for(int i = 0; i < len; i++)
  {
    h += data[i];
  }

  *(uint32_t*)out = h;
}

void sumhash32 ( const void * key, int len, uint32_t seed, void * out )
{
  uint32_t h = seed;

  const uint32_t * data = (const uint32_t*)key;

  for(int i = 0; i < len/4; i++)
  {
    h += data[i];
  }

  *(uint32_t*)out = h;
}

void DoNothingHash ( const void *, int, uint32_t, void * )
{
}

//-----------------------------------------------------------------------------
// One-byte-at-a-time hash based on Murmur's mix

uint32_t MurmurOAAT ( const void * key, int len, uint32_t seed )
{
  const uint8_t * data = (const uint8_t*)key;

  uint32_t h = seed;

  for(int i = 0; i < len; i++)
  {
    h ^= data[i];
    h *= 0x5bd1e995;
    h ^= h >> 15;
  }

  return h;
}

void MurmurOAAT_test ( const void * key, int len, uint32_t seed, void * out )
{
	*(uint32_t*)out = MurmurOAAT(key,len,seed);
}

//----------------------------------------------------------------------------

void FNV ( const void * key, int len, uint32_t seed, void * out )
{
  unsigned int h = seed;

  const uint8_t * data = (const uint8_t*)key;

  h ^= BIG_CONSTANT(2166136261);

  for(int i = 0; i < len; i++)
  {
    h ^= data[i];
    h *= 16777619;
  }

  *(uint32_t*)out = h;
}

//-----------------------------------------------------------------------------

uint32_t x17 ( const void * key, int len, uint32_t h ) 
{
  const uint8_t * data = (const uint8_t*)key;
    
  for(int i = 0; i < len; ++i) 
  {
        h = 17 * h + (data[i] - ' ');
    }

    return h ^ (h >> 16);
}

//-----------------------------------------------------------------------------

void Bernstein ( const void * key, int len, uint32_t seed, void * out ) 
{
  const uint8_t * data = (const uint8_t*)key;
    
  for(int i = 0; i < len; ++i) 
  {
        seed = 33 * seed + data[i];
    }

  *(uint32_t*)out = seed;
}

//-----------------------------------------------------------------------------
// Crap8 hash from http://www.team5150.com/~andrew/noncryptohashzoo/Crap8.html

uint32_t Crap8( const uint8_t *key, uint32_t len, uint32_t seed ) {
  #define c8fold( a, b, y, z ) { p = (uint32_t)(a) * (uint64_t)(b); y ^= (uint32_t)p; z ^= (uint32_t)(p >> 32); }
  #define c8mix( in ) { h *= m; c8fold( in, m, k, h ); }

  const uint32_t m = 0x83d2e73b, n = 0x97e1cc59, *key4 = (const uint32_t *)key;
  uint32_t h = len + seed, k = n + len;
  uint64_t p;

  while ( len >= 8 ) { c8mix(key4[0]) c8mix(key4[1]) key4 += 2; len -= 8; }
  if ( len >= 4 ) { c8mix(key4[0]) key4 += 1; len -= 4; }
  if ( len ) { c8mix( key4[0] & ( ( 1 << ( len * 8 ) ) - 1 ) ) }
  c8fold( h ^ k, n, k, k )
  return k;
}

void Crap8_test ( const void * key, int len, uint32_t seed, void * out )
{
  *(uint32_t*)out = Crap8((const uint8_t*)key,len,seed);
}
