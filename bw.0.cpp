#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <stdlib.h>
using namespace std;


//-------- global constants----------
const int N = 100;
const double G = 9.81; // m/s^2
const double DT = 0.01; // s^2 -> step size
const double K = 300.234; // N/m
const double V = 10; // m/s
const double Lx=150; // size of the box in X 
const double Ly=150;// size of the box in Y
const double long_time=100000; // number of time steps
unsigned t0, t1; // relog

//-------- body class ------------
struct Body {
  double Rxold = 0, Ryold = 0, Rx = 0, Ry = 0, Vx = 0, Vy = 0, Fx = 0, Fy = 0;
  double mass = 0;
  double rad = 0;
  void arranque(double dt);
  void timestep(double dt);
};

void Body::arranque(double dt)
{
  Rxold = Rx - dt*Vx + Fx*dt*dt/(2*mass);   
  Ryold = Ry - dt*Vy + Fy*dt*dt/(2*mass);   
}

void Body::timestep(double dt)
{
  double tmp;

  tmp = Rx;
  Rx = 2*Rx - Rxold + Fx*dt*dt/(mass);   
  Vx = (Rx - Rxold)/(2*dt);
  Rxold = tmp;

  tmp = Ry;
  Ry = 2*Ry - Ryold + Fy*dt*dt/(mass);   
  Vy = (Ry - Ryold)/(2*dt);
  Ryold = tmp;
}

// ------ function declarations-----------

void set_masses(Body bodies[]);
 
void compute_forces(Body bodies[]);
void start(Body bodies[], double dt);
void evolve(Body bodies[], double dt);

void init_gnuplot(void);
void print_to_gnuplot(Body bodies[]);

// --------- (( MAIN )) ----------------
int main(void)
{ 
  int a,b;
  srand(0); //random seed

  std::ofstream fout("datos.txt");
  
  Body bodies[N];
  for(int i=0; i < N; ++i){
   a= rand()%10;
   b= rand()%10; 
   bodies[i].rad = 1;
   bodies[i].Rx = rand()%100; // random position  in x
   bodies[i].Ry = rand()%100; // random position in y
   bodies[i].Vx = a/(sqrt(a*a+b*b))*V*pow(-1,i) ; // velocidades aleatorias 
   bodies[i].Vy = b/(sqrt(a*a+b*b))*V*pow(-1,i) ;
  }

  set_masses(bodies);
  compute_forces(bodies);

  init_gnuplot();

  start(bodies, DT);

   t0=clock(); // init clock

  for (int it = 0; it < long_time; ++it) {

    fout << DT*it << " , " << 0.5*bodies[0].mass*(bodies[0].Vx*bodies[0].Vx+bodies[0].Vy*bodies [0].Vy)+bodies[0].mass*G*bodies[0].Ry << std::endl; 
  
    compute_forces(bodies);
    evolve(bodies, DT);
    print_to_gnuplot(bodies);
  }

  fout.close();

  t1 = clock(); // finsh clock
  double time = (double(t1-t0)/CLOCKS_PER_SEC); 
  cout << "Time: " << time << endl; // print the time
  return 0;
}
//-------Mass asignation ----------

void set_masses(Body bodies[])
{
  int ii;
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].mass = 1 + double(rand())/RAND_MAX;
  }  
}


//--------- function definitions-----------
void compute_forces(Body bodies[])
{
  int ii;
  double delta;

  // ---------reset forces----------
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].Fx = 0.0;
    bodies[ii].Fy = 0.0;
  }

  
  //-------- add force with bottom wall y -----------
  for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - bodies[ii].Ry;
    if (delta > 0) {
      bodies[ii].Fy += K*delta;
     }
   }

   //------ add force with bottom wall x left------
  for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - bodies[ii].Rx;
    if (delta > 0) {
      bodies[ii].Fx += K*delta;
     }
   }
   

 
    //----- add force with bottom wall x rigth ----- 
   for (ii = 0; ii < N; ++ii) {
    delta = bodies[ii].rad - (Lx-bodies[ii].Rx);
     if (delta > 0) {
       bodies[ii].Fx += -K*delta;
      }
    }

       //----- add force with bottom wall Y up ----- 
   for (ii = 0; ii < N; ++ii) {
    delta = (bodies[ii].rad+bodies[ii].Ry)-Ly;
     if (delta > 0) {
       bodies[ii].Fy += -K*delta;
      }
    }

  // ----- fuerza with other bodies -----
  int jj;
  double Rijx, Rijy, Rij, Fx, Fy;
  for (ii = 0; ii < N; ++ii) {
    for (jj = ii+1; jj < N; ++jj) {
      Rijx = bodies[ii].Rx - bodies[jj].Rx;
      Rijy = bodies[ii].Ry - bodies[jj].Ry;
      Rij = std::sqrt(Rijx*Rijx + Rijy*Rijy);
      delta = bodies[ii].rad + bodies[jj].rad - Rij;
      if (delta > 0) {
	Fx = K*delta*Rijx/Rij;
	Fy = K*delta*Rijy/Rij;
	bodies[ii].Fx += Fx;
	bodies[ii].Fy += Fy;
	bodies[jj].Fx -= Fx;
	bodies[jj].Fy -= Fy;
      }
    }
  }  


}


void start(Body bodies[], double dt)
{
  int ii;
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].arranque(dt);
  }
}

void evolve(Body bodies[], double dt)
{
  int ii;
  for (ii = 0; ii < N; ++ii) {
    bodies[ii].timestep(dt);
  }
}
//---- Graphing ----------

void init_gnuplot(void)
{
  std::cout << "set size ratio -1" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange [0:1]" << std::endl;
   std::cout << "set xrange [0:" << Lx << "]" << std::endl; 
  std::cout  << "set yrange [0: "<< Ly << "]" << std::endl;
  std::cout << "unset key" << std::endl; //---> quita nombre de los  graficos.
}


void print_to_gnuplot(Body bodies[])
{
  std::cout << "plot "; 
  for (int ii = 0; ii < N; ++ii) {
    std::cout << bodies[ii].Rx << " + " << bodies[ii].rad << "*cos(t*2*pi) , "
	      << bodies[ii].Ry << " + " << bodies[ii].rad << "*sin(t*2*pi) , ";
  }
  std::cout << " 0, 0"; 
  std::cout << std::endl;
}
