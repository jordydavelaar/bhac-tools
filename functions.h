/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 * A list of all RAPTOR functions.
 */

#include "definitions.h"
#include "model_definitions.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H


// UNIFORM.C
////////////

void uniform_data_array(double TIME_INIT);


// CORE.C
/////////

//////////////////////

// See grmonty paper by Dolence et al.
// HARM model internal utilities
void init_model();

void set_units(double);

void init_grmhd_data(char *fname);

void init_storage();

void Xtoijk(double *X, int *i, int *j, int *k, double *del);

// void get_fluid_params(double X[4], double *Ne, double *Thetae, double *B,
//                      double *B_u, double Ucon[4], int *IN_VOLUME);
int get_fluid_params(double X[NDIM], struct GRMHD *modvar);
// IO


// GRMATH.C
///////////

void get_spatial_lc(double lc[NDIM][NDIM][NDIM]);

void get_full_lc(double lc[NDIM][NDIM][NDIM]);


double get_r(double X_u[4]);

// Lowers the index of the contravariant vector V_u, storing the results in
// a covariant one (V_d), based on the metric at position X_u
void lower_index(double X_u[4], double V_u[4], double V_d[4]);

// Lowers two indices of a rank (2, 0) tensor
void lower_two_indices(double N_uu[4][4], double N_dd[4][4], double X_u[4]);

// Lowers the index of a contravariant vector V_u in BL coordinates.
void BL_lower_index(double X_u[4], double V_u[4], double V_d[4]);

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
void raise_index(double X_u[4], double V_d[4], double V_u[4]);

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
// Needed for CKS coordinates
void raise_index_KS(double X_u[4], double V_d[4], double V_u[4]);

// Adjusts y[4] = U_u[0] so that y describes a lightray/null geodesic
void normalize_null(double X_u[4], double U_u[4]);

// Returns the norm of U_u, which is the scalar g_dd[a][b] * U_u[a] * U_u[b]
double four_velocity_norm(double X_u[4], double U_u[4]);

// Returns the inner product of vectors A and B, i.e. A_u B_d
double inner_product(double *X_u, double *A_u, double *B_u);

// Transform a contravariant vector from BL to KS coordinates
void BL_to_KS_u(double *BLphoton_u, double *KSphoton_u);

// Transform a contravariant vector from KS to BL coordinates
void KS_to_BL_u(double *KSphoton_u, double *BLphoton_u);

// Convert KS to CKS coordinates
void KS_to_CKS(double *X_KS_u, double *X_CKS_u);

// Convert CKS to KS coordinates
void CKS_to_KS(double *X_CKS_u, double *X_KS_u);

// Transform a contravariant vector from KS to CKS coordinates
void KS_to_CKS_u(double *KScoords, double *CKScoords);

// Return the photon frequency in the co-moving frame of the plasma
double freq_in_plasma_frame(double Uplasma_u[4], double k_d[4]);

// Angle between k_u and B_u in the plasma frame
double pitch_angle(double *X_u, double *k_u, double *B_u, double *Uplasma_u);

// void f_tetrad_to_stokes(double Iinv, double Iinv_pol, double complex
// f_tetrad_u[], double complex S_A[4]);

// void stokes_to_f_tetrad(double complex S_A[], double *Iinv, double *Iinv_pol,
// double complex f_tetrad_u[4]);

// void construct_U_vector( double X_u[], double U_u[]);

// METRIC.C
///////////

// Computes the metric at location X
void metric_dd(double X_u[4], double g_dd[4][4]);

// Computes the inverse metric at location X
void metric_uu(double X_u[4], double g_uu[4][4]);

// Computes the inverse metric at location X
void metric_KS_uu(double X_u[4], double g_uu[4][4]);

// Computes the Christoffel symbols at location X numerically (general metric)
void connection_num_udd(double X_u[4], double gamma_udd[4][4][4]);

// Computes the Christoffel symbols at location X based on an exact metric
void connection_udd(double X_u[4], double gamma_udd[4][4][4]);

// This function initializes a single 'superphoton' or light ray.
void initialize_photon(double alpha, double beta, double k_u[4], double t_init);

// Transformation functions
double Xg2_approx_rand(double Xr2);

double Ug2_approx_rand(double Ur2, double Xg2);


#endif // FUNCTIONS_H
