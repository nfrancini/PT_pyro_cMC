#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    double Jzz, Jpm, Jpmpm, Jzpm;
    double j1, j2, j3, j4;
    double j, j_ising, j_pd, j_dm;

    // Prompt the user to enter four real numbers
    printf("Enter Jzz: ");
    scanf("%lf", &Jzz);

    printf("Enter J+-: ");
    scanf("%lf", &Jpm);

    printf("Enter J+-+-: ");
    scanf("%lf", &Jpmpm);

    printf("Enter Jz+-: ");
    scanf("%lf", &Jzpm);

    // Display the entered numbers
    printf("\nYou entered the following numbers:\n");
    printf("Jzz:\t%.13lf\n", Jzz);
    printf("J+-:\t%.13lf\n", Jpm);
    printf("J+-+-:\t%.13lf\n", Jpmpm);
    printf("Jz+-:\t%.13lf\n", Jzpm);

    // First conversion from {Jzz, J+-, J+-, Jz+-} to {J1,J2,J3,J4}
    j1=(-Jzz+4*Jpm+2*Jpmpm+2*sqrt(2)*Jzpm)/3.0;
    j2=(Jzz-4*Jpm+4*Jpmpm+4*sqrt(2)*Jzpm)/3.0;
    j3=(-Jzz-2*Jpm-4*Jpmpm+2*sqrt(2)*Jzpm)/3.0;
    j4=(-Jzz-2*Jpm+2*Jpmpm-sqrt(2)*Jzpm)/3.0;

    // Display first conversion number
    printf("\nFirst conversion:\n");
    printf("J1:\t%.13lf\n", j1);
    printf("J2:\t%.13lf\n", j2);
    printf("J3:\t%.13lf\n", j3);
    printf("J4:\t%.13lf\n", j4);

    // Direct converion from {Jzz, J+-, J+-, Jz+-} to {J,J_ising,J_PD,J_DM}
    j=2*(Jpm + Jpmpm);
    j_ising=Jzz -10*Jpm -2*Jpmpm +4*sqrt(2)*Jzpm;
    j_pd=-8*Jpm -4*Jpmpm +4*sqrt(2)*Jzpm;
    j_dm=-4*Jpm + sqrt(2)*Jzpm;

    // Display direct conversion number
    printf("\nDirect conversion:\n");
    printf("J:\t%.13lf\n", j);
    printf("J_isi:\t%.13lf\n", j_ising);
    printf("J_PD:\t%.13lf\n", j_pd);
    printf("J_DM:\t%.13lf\n", j_dm);

    return 0;
}
