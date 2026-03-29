#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "polyvec_accel.h"

int check_result(const int16_t *hw_res, const int16_t *sw_res, int size, const char* op_name) {
    int errors = 0;
    for(int i = 0; i < size; i++) {
        if(hw_res[i] != sw_res[i]) {
            if (errors < 10) { 
                printf("  ERROR [%s]: idx[%d] HW = %d, SW (Expected) = %d\n", op_name, i, hw_res[i], sw_res[i]);
            }
            errors++;
        }
    }
    if (errors > 0) printf("  -> Total of errors in %s: %d\n", op_name, errors);
    else printf("  -> %s: Passed!\n", op_name);
    return errors;
}

int main() {
    int total_errors = 0;

    int16_t *mem_pv_r = (int16_t *)malloc(KYBER_N * KYBER_K * sizeof(int16_t));
    int16_t *mem_pv_a = (int16_t *)malloc(KYBER_N * KYBER_K * sizeof(int16_t));
    int16_t *mem_pv_b = (int16_t *)malloc(KYBER_N * KYBER_K * sizeof(int16_t));
    
    int16_t *mem_p_r  = (int16_t *)malloc(KYBER_N * sizeof(int16_t));
    int16_t *mem_p_a  = (int16_t *)malloc(KYBER_N * sizeof(int16_t));
    int16_t *mem_p_b  = (int16_t *)malloc(KYBER_N * sizeof(int16_t));

    polyvec sw_pv_r, sw_pv_a, sw_pv_b;
    poly sw_p_r, sw_p_a, sw_p_b;

    for(int i = 0; i < KYBER_K; i++) {
        for(int j = 0; j < KYBER_N; j++) {
            sw_pv_a.vec[i].coeffs[j] = rand() % KYBER_Q;
            sw_pv_b.vec[i].coeffs[j] = rand() % KYBER_Q;
        }
    }
    for(int i = 0; i < KYBER_N; i++) {
        sw_p_a.coeffs[i] = rand() % KYBER_Q;
        sw_p_b.coeffs[i] = rand() % KYBER_Q;
    }

    // ==========================================
    // POLYVEC TESTS
    // ==========================================

    // OP 0: POLYVEC_NTT
    printf("OP 0: POLYVEC_NTT\n");
    sw_pv_r = sw_pv_a; 
    memcpy(mem_pv_r, &sw_pv_r, KYBER_N * KYBER_K * sizeof(int16_t));
    polyvec_ntt(&sw_pv_r);
    polyvec_accel(OP_POLYVEC_NTT, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b); 
    total_errors += check_result(mem_pv_r, (int16_t*)&sw_pv_r, KYBER_N * KYBER_K, "OP_POLYVEC_NTT");

    // OP 1: POLYVEC_INVNTT_TOMONT
    printf("OP 1: POLYVEC_INVNTT_TOMONT\n");
    sw_pv_r = sw_pv_a; 
    memcpy(mem_pv_r, &sw_pv_r, KYBER_N * KYBER_K * sizeof(int16_t));
    polyvec_invntt_tomont(&sw_pv_r); 
    polyvec_accel(OP_POLYVEC_INVNTT_TOMONT, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b);
    total_errors += check_result(mem_pv_r, (int16_t*)&sw_pv_r, KYBER_N * KYBER_K, "OP_POLYVEC_INVNTT_TOMONT");

    // OP 2: POLYVEC_BASEMUL_ACC_MONTGOMERY
    printf("OP 2: POLYVEC_BASEMUL_ACC_MONTGOMERY\n");
    memcpy(mem_pv_a, &sw_pv_a, KYBER_N * KYBER_K * sizeof(int16_t));
    memcpy(mem_pv_b, &sw_pv_b, KYBER_N * KYBER_K * sizeof(int16_t));
    polyvec_basemul_acc_montgomery(&sw_p_r, &sw_pv_a, &sw_pv_b); 
    polyvec_accel(OP_POLYVEC_BASEMUL_ACC_MONTGOMERY, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b); 
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLYVEC_BASEMUL_ACC_MONTGOMERY");

    // OP 3: OP_POLYVEC_ADD
    printf("OP 3: OP_POLYVEC_ADD\n");
    memcpy(mem_pv_a, &sw_pv_a, KYBER_N * KYBER_K * sizeof(int16_t));
    memcpy(mem_pv_b, &sw_pv_b, KYBER_N * KYBER_K * sizeof(int16_t));
    polyvec_add(&sw_pv_r, &sw_pv_a, &sw_pv_b); 
    polyvec_accel(OP_POLYVEC_ADD, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b);
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLYVEC_ADD");

    // OP 4: POLYVEC_REDUCE
    printf("OP 4: POLYVEC_REDUCE\n");
    sw_pv_r = sw_pv_a; 
    memcpy(mem_pv_r, &sw_pv_r, KYBER_N * KYBER_K * sizeof(int16_t));
    polyvec_reduce(&sw_pv_r); 
    polyvec_accel(OP_POLYVEC_REDUCE, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b);
    total_errors += check_result(mem_pv_r, (int16_t*)&sw_pv_r, KYBER_N * KYBER_K, "OP_POLYVEC_REDUCE");

    // ==========================================
    // POLY TESTS
    // ==========================================

    // OP 5: POLY_TOMONT
    printf("OP 5: POLY_TOMONT\n");
    sw_p_r = sw_p_a; 
    memcpy(mem_p_r, &sw_p_r, KYBER_N * sizeof(int16_t));
    poly_tomont(&sw_p_r); 
    polyvec_accel(OP_POLY_TOMONT, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b); 
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLY_TOMONT");

    // OP 6: POLY_INVNTT_TOMONT
    printf("OP 6: POLY_INVNTT_TOMONT\n");
    sw_p_r = sw_p_a; 
    memcpy(mem_p_r, &sw_p_r, KYBER_N * sizeof(int16_t));
    poly_invntt_tomont(&sw_p_r); 
    polyvec_accel(OP_POLY_INVNTT_TOMONT, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b); 
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLY_INVNTT_TOMONT");

    // OP 7: OP_POLY_ADD
    printf("OP 7: OP_POLY_ADD\n");

    memcpy(mem_p_a, &sw_p_a, KYBER_N * sizeof(int16_t));
    memcpy(mem_p_b, &sw_p_b, KYBER_N * sizeof(int16_t));
    poly_add(&sw_p_r, &sw_p_a, &sw_p_b); 
    polyvec_accel(OP_POLY_ADD, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b);  
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLY_ADD");

    // OP 8: POLY_SUB
    printf("OP 8: POLY_SUB\n");
    memcpy(mem_p_a, &sw_p_a, KYBER_N * sizeof(int16_t));
    memcpy(mem_p_b, &sw_p_b, KYBER_N * sizeof(int16_t));
    poly_sub(&sw_p_r, &sw_p_a, &sw_p_b); 
    polyvec_accel(OP_POLY_SUB, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b); 
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLY_SUB");

    // OP 9: POLY_REDUCE
    printf("OP 9: POLY_REDUCE\n");
    sw_p_r = sw_p_a; 
    memcpy(mem_p_r, &sw_p_r, KYBER_N * sizeof(int16_t));
    poly_reduce(&sw_p_r); 
    polyvec_accel(OP_POLY_REDUCE, mem_pv_r, mem_pv_a, mem_pv_b, mem_p_r, mem_p_a, mem_p_b); 
    total_errors += check_result(mem_p_r, (int16_t*)&sw_p_r, KYBER_N, "OP_POLY_REDUCE");

    free(mem_pv_r); free(mem_pv_a); free(mem_pv_b);
    free(mem_p_r);  free(mem_p_a);  free(mem_p_b);

    if(total_errors == 0) {
        printf("Success.\n");
        return 0; 
    } else {
        printf("Failed: Found %d errors in total.\n", total_errors);
        return 1;
    }
}