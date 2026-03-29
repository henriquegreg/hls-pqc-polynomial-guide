#include <stdint.h>
#include <string.h>

#define OP_POLYVEC_NTT 0
#define OP_POLYVEC_INVNTT_TOMONT 1
#define OP_POLYVEC_BASEMUL_ACC_MONTGOMERY 2
#define OP_POLYVEC_ADD 3 
#define OP_POLYVEC_REDUCE 4

#define OP_POLY_TOMONT 5
#define OP_POLY_INVNTT_TOMONT 6
#define OP_POLY_ADD 7 
#define OP_POLY_SUB 8 
#define OP_POLY_REDUCE 9

void polyvec_accel(
    uint8_t op,
    int16_t *mem_pv_r, const int16_t *mem_pv_a, const int16_t *mem_pv_b,
    int16_t *mem_p_r,  const int16_t *mem_p_a,  const int16_t *mem_p_b
);
