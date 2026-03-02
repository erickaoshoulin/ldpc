#include "ldpc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *name;
    int n;
    int m;
    int dv;
} matrix_profile_t;

static int write_performance_plots(void) {
    int rc1 = system("python3 scripts/plot_results.py >/dev/null 2>&1");
    int rc2 = system("python3 scripts/plot_matrix_results.py >/dev/null 2>&1");
    if (rc1 != 0 || rc2 != 0) {
        fprintf(stderr, "warning: failed to generate one or more figures (requires python3)\n");
        return -1;
    }
    return 0;
}

static void run_eval(const ldpc_matrix_t *h, ldpc_algorithm_t alg, const ldpc_decoder_params_t *p,
                     float snr_db, int frames, uint32_t seed,
                     double *ber, double *fer, double *avg_iters) {
    uint8_t *tx = (uint8_t *)calloc((size_t)h->n, sizeof(uint8_t));
    uint8_t *rx = (uint8_t *)calloc((size_t)h->n, sizeof(uint8_t));
    float *llr = (float *)calloc((size_t)h->n, sizeof(float));
    if (!tx || !rx || !llr) {
        fprintf(stderr, "allocation failed\n");
        exit(1);
    }

    long long bit_err = 0;
    long long frame_err = 0;
    long long iter_sum = 0;

    const float rate = (float)(h->n - h->m) / (float)h->n;
    for (int f = 0; f < frames; ++f) {
        memset(tx, 0, (size_t)h->n);
        ldpc_awgn_bpsk_llr(tx, h->n, snr_db, rate, &seed, llr);
        ldpc_decode_result_t r = ldpc_decode(h, llr, alg, p, rx);

        int fe = 0;
        for (int i = 0; i < h->n; ++i) {
            bit_err += rx[i] != tx[i];
            fe |= rx[i] != tx[i];
        }
        frame_err += fe;
        iter_sum += r.iterations_used;
    }

    *ber = (double)bit_err / (double)(frames * h->n);
    *fer = (double)frame_err / (double)frames;
    *avg_iters = (double)iter_sum / (double)frames;

    free(tx);
    free(rx);
    free(llr);
}

int main(int argc, char **argv) {
    int n = 512;
    int m = 64;
    int dv = 6;
    int max_iters = 20;
    int frames = 300;
    int group_size = 8;
    float alpha = 0.80f;
    float beta = 0.15f;
    float damping = 0.20f;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) m = atoi(argv[2]);
    if (argc > 3) dv = atoi(argv[3]);
    if (argc > 4) frames = atoi(argv[4]);

    const matrix_profile_t profiles[] = {
        {"custom", n, m, dv},
        {"ieee_80211n_r12", 648, 324, 6},
        {"ieee_80211n_r23", 648, 216, 6},
        {"ieee_80211n_r34", 648, 162, 6},
    };

    ldpc_decoder_params_t conv = {
        .max_iters = max_iters,
        .alpha = alpha,
        .alpha_final = alpha,
        .beta = beta,
        .group_size = 1,
        .damping = 0.0f,
    };
    ldpc_decoder_params_t rmas1 = {
        .max_iters = max_iters,
        .alpha = alpha,
        .alpha_final = alpha,
        .beta = beta,
        .group_size = group_size,
        .damping = damping,
    };
    ldpc_decoder_params_t rmas2 = {
        .max_iters = max_iters,
        .alpha = alpha,
        .alpha_final = alpha,
        .beta = beta * 0.8f,
        .group_size = group_size,
        .damping = damping * 0.5f,
    };
    ldpc_decoder_params_t as = {
        .max_iters = max_iters,
        .alpha = 0.70f,
        .alpha_final = 0.95f,
        .beta = beta * 0.5f,
        .group_size = 1,
        .damping = 0.0f,
    };

    FILE *fp = fopen("results/performance.csv", "w");
    FILE *fp_matrix = fopen("results/matrix_performance.csv", "w");
    if (!fp || !fp_matrix) {
        perror("fopen");
        if (fp) fclose(fp);
        if (fp_matrix) fclose(fp_matrix);
        return 1;
    }
    fprintf(fp, "snr_db,algorithm,ber,fer,avg_iterations\n");
    fprintf(fp_matrix, "profile,snr_db,algorithm,ber,fer,avg_iterations\n");

    for (size_t pi = 0; pi < sizeof(profiles) / sizeof(profiles[0]); ++pi) {
        const matrix_profile_t *mp = &profiles[pi];
        ldpc_matrix_t h;
        if (ldpc_generate_regular(&h, mp->n, mp->m, mp->dv, 0x12345678u + (uint32_t)pi) != 0) {
            fprintf(stderr, "failed to generate H matrix for profile %s\n", mp->name);
            fclose(fp);
            fclose(fp_matrix);
            return 1;
        }

        printf("\nLDPC eval profile=%s n=%d m=%d dv=%d rate=%.4f frames=%d\n", mp->name, h.n, h.m, mp->dv,
               (double)(h.n - h.m) / h.n, frames);

        for (float snr = 4.5f; snr <= 7.0f; snr += 0.5f) {
            double ber, fer, it;
            run_eval(&h, LDPC_ALG_CONVENTIONAL, &conv, snr, frames, 0xABCDEF01u, &ber, &fer, &it);
            if (pi == 0) fprintf(fp, "%.2f,conventional,%.8e,%.8e,%.4f\n", snr, ber, fer, it);
            fprintf(fp_matrix, "%s,%.2f,conventional,%.8e,%.8e,%.4f\n", mp->name, snr, ber, fer, it);

            run_eval(&h, LDPC_ALG_RMAS1, &rmas1, snr, frames, 0x10203040u, &ber, &fer, &it);
            if (pi == 0) fprintf(fp, "%.2f,rmas1,%.8e,%.8e,%.4f\n", snr, ber, fer, it);
            fprintf(fp_matrix, "%s,%.2f,rmas1,%.8e,%.8e,%.4f\n", mp->name, snr, ber, fer, it);

            run_eval(&h, LDPC_ALG_RMAS2, &rmas2, snr, frames, 0x55667788u, &ber, &fer, &it);
            if (pi == 0) fprintf(fp, "%.2f,rmas2,%.8e,%.8e,%.4f\n", snr, ber, fer, it);
            fprintf(fp_matrix, "%s,%.2f,rmas2,%.8e,%.8e,%.4f\n", mp->name, snr, ber, fer, it);

            run_eval(&h, LDPC_ALG_AS, &as, snr, frames, 0x77AABBCCu, &ber, &fer, &it);
            if (pi == 0) fprintf(fp, "%.2f,as,%.8e,%.8e,%.4f\n", snr, ber, fer, it);
            fprintf(fp_matrix, "%s,%.2f,as,%.8e,%.8e,%.4f\n", mp->name, snr, ber, fer, it);
        }
        ldpc_free_matrix(&h);
    }

    fclose(fp);
    fclose(fp_matrix);
    write_performance_plots();

    printf("results written to results/performance.csv, results/matrix_performance.csv and figures\n");
    return 0;
}
