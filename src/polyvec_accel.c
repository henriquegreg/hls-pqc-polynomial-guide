#include "polyvec.h"
#include "poly.h"
#include "polyvec_accel.h"

void polyvec_accel(
    uint8_t op,
    int16_t *mem_pv_r, const int16_t *mem_pv_a, const int16_t *mem_pv_b,
    int16_t *mem_p_r,  const int16_t *mem_p_a,  const int16_t *mem_p_b
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=control
    #pragma HLS INTERFACE s_axilite port=op bundle=control

    #pragma HLS INTERFACE m_axi port=mem_pv_r depth=KYBER_N*KYBER_K bundle=gmem0 offset=slave \
        num_write_outstanding=16 max_write_burst_length=256
    #pragma HLS INTERFACE s_axilite port=mem_pv_r bundle=control

    #pragma HLS INTERFACE m_axi port=mem_pv_a depth=KYBER_N*KYBER_K bundle=gmem1 offset=slave \
        num_read_outstanding=16 max_read_burst_length=256
    #pragma HLS INTERFACE s_axilite port=mem_pv_a bundle=control

    #pragma HLS INTERFACE m_axi port=mem_pv_b depth=KYBER_N*KYBER_K bundle=gmem2 offset=slave \
        num_read_outstanding=16 max_read_burst_length=256
    #pragma HLS INTERFACE s_axilite port=mem_pv_b bundle=control

    #pragma HLS INTERFACE m_axi port=mem_p_r depth=KYBER_N bundle=gmem0 offset=slave \
        num_write_outstanding=16 max_write_burst_length=256
    #pragma HLS INTERFACE s_axilite port=mem_p_r bundle=control

    #pragma HLS INTERFACE m_axi port=mem_p_a depth=KYBER_N bundle=gmem1 offset=slave \
        num_read_outstanding=16 max_read_burst_length=256
    #pragma HLS INTERFACE s_axilite port=mem_p_a bundle=control

    #pragma HLS INTERFACE m_axi port=mem_p_b depth=KYBER_N bundle=gmem2 offset=slave \
        num_read_outstanding=16 max_read_burst_length=256
    #pragma HLS INTERFACE s_axilite port=mem_p_b bundle=control

    polyvec local_pv_r, local_pv_a, local_pv_b;
    poly    local_p_r,  local_p_a,  local_p_b;

    switch(op) {
        case OP_POLYVEC_NTT:
            memcpy(&local_pv_r, mem_pv_r, KYBER_N*KYBER_K * sizeof(int16_t));
            polyvec_ntt(&local_pv_r); 
            memcpy(mem_pv_r, &local_pv_r, KYBER_N*KYBER_K * sizeof(int16_t)); 
            break;
            
        case OP_POLYVEC_INVNTT_TOMONT:
            memcpy(&local_pv_r, mem_pv_r, KYBER_N*KYBER_K * sizeof(int16_t));
            polyvec_invntt_tomont(&local_pv_r);
            memcpy(mem_pv_r, &local_pv_r, KYBER_N*KYBER_K * sizeof(int16_t));
            break;
            
        case OP_POLYVEC_BASEMUL_ACC_MONTGOMERY:
            memcpy(&local_pv_a, mem_pv_a, KYBER_N*KYBER_K * sizeof(int16_t));
            memcpy(&local_pv_b, mem_pv_b, KYBER_N*KYBER_K * sizeof(int16_t));
            polyvec_basemul_acc_montgomery(&local_p_r, &local_pv_a, &local_pv_b);
            memcpy(mem_p_r, &local_p_r, KYBER_N * sizeof(int16_t));
            break;
            
        case OP_POLYVEC_ADD:
            memcpy(&local_pv_a, mem_pv_a, KYBER_N*KYBER_K * sizeof(int16_t));
            memcpy(&local_pv_b, mem_pv_b, KYBER_N*KYBER_K * sizeof(int16_t));
            polyvec_add(&local_pv_r, &local_pv_a, &local_pv_b);
            memcpy(mem_pv_r, &local_pv_r, KYBER_N*KYBER_K * sizeof(int16_t));
            break;
            
        case OP_POLYVEC_REDUCE:
            memcpy(&local_pv_r, mem_pv_r, KYBER_N*KYBER_K * sizeof(int16_t));
            polyvec_reduce(&local_pv_r);
            memcpy(mem_pv_r, &local_pv_r, KYBER_N*KYBER_K * sizeof(int16_t));
            break;

        case OP_POLY_TOMONT:
            memcpy(&local_p_r, mem_p_r, KYBER_N * sizeof(int16_t));
            poly_tomont(&local_p_r);
            memcpy(mem_p_r, &local_p_r, KYBER_N * sizeof(int16_t));
            break;
            
        case OP_POLY_INVNTT_TOMONT:
        {
            memcpy(&local_p_r, mem_p_r, KYBER_N * sizeof(int16_t));
            
            poly tmp = local_p_r;
            #pragma HLS ARRAY_PARTITION variable=tmp.coeffs type=cyclic factor=16 dim=1
            
            poly_invntt_tomont(&tmp);
            
            local_p_r = tmp;
            
            memcpy(mem_p_r, &local_p_r, KYBER_N * sizeof(int16_t));
            break;
        }
            
        case OP_POLY_ADD:
            memcpy(&local_p_a, mem_p_a, KYBER_N * sizeof(int16_t));
            memcpy(&local_p_b, mem_p_b, KYBER_N * sizeof(int16_t));
            poly_add(&local_p_r, &local_p_a, &local_p_b);
            memcpy(mem_p_r, &local_p_r, KYBER_N * sizeof(int16_t));
            break;
            
        case OP_POLY_SUB:
            memcpy(&local_p_a, mem_p_a, KYBER_N * sizeof(int16_t));
            memcpy(&local_p_b, mem_p_b, KYBER_N * sizeof(int16_t));
            poly_sub(&local_p_r, &local_p_a, &local_p_b);
            memcpy(mem_p_r, &local_p_r, KYBER_N * sizeof(int16_t));
            break;
            
        case OP_POLY_REDUCE:
            memcpy(&local_p_r, mem_p_r, KYBER_N * sizeof(int16_t));
            poly_reduce(&local_p_r);
            memcpy(mem_p_r, &local_p_r, KYBER_N * sizeof(int16_t));
            break;

        default:
            break;
    }
}