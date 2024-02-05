
#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

char GRMHD_FILE[256];

double MBH, M_UNIT, TIME_INIT, INCLINATION;
double R_HIGH, R_LOW;
double FREQS_PER_DEC, FREQ_MIN, FREQ_MAX;

double SOURCE_DIST; // Distance to M87 (cm); for Sgr A* use (2.47e22)

int IMG_WIDTH, IMG_HEIGHT;
double CAM_SIZE_X, CAM_SIZE_Y;
double STEPSIZE;

int main(int argc, char *argv[]) {

    // INPUT FILE
    /////////////

    fprintf(stderr, "\nInitializing...\n");
    // INITIALIZE MODEL
    ///////////////////
    sscanf(argv[1], "%s", GRMHD_FILE);
    sscanf(argv[2], "%lf", &TIME_INIT);
    fprintf(stderr,"t=%e\n",TIME_INIT);
    // Initialize HARM2D grmhd model
    // Note: this sets the black hole spin 'a'
    init_model();

    uniform_data_array(TIME_INIT);

    return 0;
}
