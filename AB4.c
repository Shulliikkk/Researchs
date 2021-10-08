#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Density.h"
#define Ro(r) interpol(h, ro, l, (r))

double interpol(double* arr_x, double* arr_y, int len, double x){
  if (x >= arr_x[0] && x <= arr_x[len - 1]){
    double x0, x1, y0, y1;
    for (int i = 0; i < len - 1; i++){
      if (x >= arr_x[i] && x <= arr_x[i + 1]) {
        x0 = arr_x[i]; x1 = arr_x[i + 1];
        y0 = arr_y[i]; y1 = arr_y[i + 1];
      }
    }
    return((y1 - y0)/(x1 - x0)*(x - x0) + y0);
  }
  else {
    printf("Out of range");
    exit(1);
  }
}

int main() {
  FILE* f = fopen("y(x).txt", "w");
  FILE* f2 = fopen("h(t).txt", "w");

  double dt; printf("dt = "); scanf("%lf", &dt);
  double T; printf("T = "); scanf("%lf", &T);

  double M = 5.9722E24;
  double R = 6371E3;
  double G = 6.67430151515E-11;
  double mu = - G * M;

  double phi0; printf("phi0 = "); scanf("%lf", &phi0);
  double h0; printf("h0 = "); scanf("%lf", &h0);
  double vx0; printf("vx0 = "); scanf("%lf", &vx0);
  double vy0; printf("xy0 = "); scanf("%lf", &vy0); //sqrt(-mu/(R+h0))


  double alph; printf("alph(S/m) = "); scanf("%lf", &alph);
  double Cx = 2.5;
  const double k = Cx*alph/2 ;

  double x[5], y[5];
  double vx[5], vy[5];
  double ax[5], ay[5];
  double r;

  double t = 0;

  x[0] = (R + h0) * cos(phi0);
  y[0] = (R + h0) * sin(phi0);
  r = x[0]*x[0] + y[0]*y[0];
  vx[0] = vx0;
  vy[0] = vy0;
  ax[0] = mu/r * (x[0]/sqrt(r)) - k * Ro(sqrt(r) - R) * vx[0]*sqrt(vx[0]*vx[0] + vy[0]*vy[0]);
  ay[0] = mu/r * (y[0]/sqrt(r)) - k * Ro(sqrt(r) - R) * vy[0]*sqrt(vx[0]*vx[0] + vy[0]*vy[0]);
  fprintf(f, "%lf %lf\n", x[0], y[0]);
  t += dt/3600;
  fprintf(f2, "%lf %lf\n", t, sqrt(r) - R);

  x[1] = x[0] + vx[0]*dt;
  y[1] = y[0] + vy[0]*dt;
  r = x[1]*x[1] + y[1]*y[1];
  vx[1] = vx[0] + ax[0]*dt;
  vy[1] = vy[0] + ay[0]*dt;
  ax[1] = mu/r * (x[1]/sqrt(r)) - k * Ro(sqrt(r) - R) * vx[1]*sqrt(vx[1]*vx[1] + vy[1]*vy[1]);
  ay[1] = mu/r * (y[1]/sqrt(r)) - k * Ro(sqrt(r) - R) * vy[1]*sqrt(vx[1]*vx[1] + vy[1]*vy[1]);
  fprintf(f, "%lf %lf\n", x[1], y[1]);
  t += dt/3600;
  fprintf(f2, "%lf %lf\n", t, sqrt(r) - R);

  x[2] = x[1] + (3./2 * vx[1] - 1./2 * vx[0])*dt;
  y[2] = y[1] + (3./2 * vy[1] - 1./2 * vy[0])*dt;
  r = x[2]*x[2] + y[2]*y[2];
  vx[2] = vx[1] + (3./2 * ax[1] - 1./2 * ax[0])*dt;
  vy[2] = vy[1] + (3./2 * ay[1] - 1./2 * ay[0])*dt;
  ax[2] = mu/r * (x[2]/sqrt(r)) - k * Ro(sqrt(r) - R) * vx[2]*sqrt(vx[2]*vx[2] + vy[2]*vy[2]);
  ay[2] = mu/r * (y[2]/sqrt(r)) - k * Ro(sqrt(r) - R) * vy[2]*sqrt(vx[2]*vx[2] + vy[2]*vy[2]);
  fprintf(f, "%lf %lf\n", x[2], y[2]);
  t += dt/3600;
  fprintf(f2, "%lf %lf\n", t, sqrt(r) - R);

  x[3] = x[2] + (23./12 * vx[2] - 16./12 * vx[1] + 5./12 * vx[0])*dt;
  y[3] = y[2] + (23./12 * vy[2] - 16./12 * vy[1] + 5./12 * vy[0])*dt;
  r = x[3]*x[3] + y[3]*y[3];
  vx[3] = vx[2] + (23./12 * ax[2] - 16./12 * ax[1] + 5./12 * ax[0])*dt;
  vy[3] = vy[2] + (23./12 * ay[2] - 16./12 * ay[1] + 5./12 * ay[0])*dt;
  ax[3] = mu/r * (x[3]/sqrt(r)) - k * Ro(sqrt(r) - R) * vx[3]*sqrt(vx[3]*vx[3] + vy[3]*vy[3]);
  ay[3] = mu/r * (y[3]/sqrt(r)) - k * Ro(sqrt(r) - R) * vy[3]*sqrt(vx[3]*vx[3] + vy[3]*vy[3]);
  fprintf(f, "%lf %lf\n", x[3], y[3]);
  t += dt/3600;
  fprintf(f2, "%lf %lf\n", t, sqrt(r) - R);

  for (int step = 1; step <= 24*3600*T/dt; step ++) {
    x[4] = x[3] + (55./24 * vx[3] - 59./24 * vx[2] + 37./24 * vx[1] - 9./24 * vx[0])*dt;
    y[4] = y[3] + (55./24 * vy[3] - 59./24 * vy[2] + 37./24 * vy[1] - 9./24 * vy[0])*dt;
    r = x[4]*x[4] + y[4]*y[4];
    if (sqrt(r) - R <= 0) {
      break;
    }
    vx[4] = vx[3] + (55./24 * ax[3] - 59./24 * ax[2] + 37./24 * ax[1] - 9./24 * ax[0])*dt;
    vy[4] = vy[3] + (55./24 * ay[3] - 59./24 * ay[2] + 37./24 * ay[1] - 9./24 * ay[0])*dt;
    ax[4] = mu/r * (x[4]/sqrt(r)) - k * Ro(sqrt(r) - R) * vx[4]*sqrt(vx[4]*vx[4] + vy[4]*vy[4]);
    ay[4] = mu/r * (y[4]/sqrt(r)) - k * Ro(sqrt(r) - R) * vy[4]*sqrt(vx[4]*vx[4] + vy[4]*vy[4]);
    fprintf(f, "%lf %lf\n", x[4], y[4]);
    t += dt/3600;
    fprintf(f2, "%lf %lf\n", t, sqrt(r) - R);

    int i;
    for (i = 0; i < 5; i++){
      x[i] = x[i + 1];
      y[i] = y[i + 1];
      vx[i] = vx[i + 1];
      vy[i] = vy[i + 1];
      ax[i] = ax[i + 1];
      ay[i] = ay[i + 1];
    }
  
  }
  system("python plotter.py");
}
