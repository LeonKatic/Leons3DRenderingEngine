#include <SDL2/SDL.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_scancode.h>
#include <SDL2/SDL_stdinc.h>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

#include <vector>
using namespace std;

class Tocka {
public:
  int id;
  double x;
  double y;
  double z;
  Tocka() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
  }
  Tocka(double xkor, double ykor) {
    this->x = xkor;
    this->y = ykor;
    this->z = 0;
  }
  Tocka(double xkor, double ykor, double zkor) {
    this->x = xkor;
    this->y = ykor;
    this->z = zkor;
  }
  void printTocka() {
    cout << this->x << " " << this->y << " " << this->z << endl;
  }
};

class Vektor {
public:
  double x;
  double y;
  double z;
  Vektor() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
  }
  Vektor(double xkor, double ykor, double zkor) {
    this->x = xkor;
    this->y = ykor;
    this->z = zkor;
  }

  Vektor operator+(const Vektor &other) const {
    return Vektor(this->x + other.x, this->y + other.y, this->z + other.z);
  }
  Vektor operator-(const Vektor &other) const {
    return Vektor(this->x - other.x, this->y - other.y, this->z - other.z);
  }

  double operator*(const Vektor &other) const {
    return this->x * other.x + this->y * other.y + this->z * other.z;
  }

  Vektor operator*(double a) const {
    return Vektor(a * this->x, a * this->y, a * this->z);
  }

  void printVektor() {
    cout << this->x << " " << this->y << " " << this->z << endl;
  }
};
Vektor operator*(double a, const Vektor &v) {
  return v * a; // reuse member operator*
}

// cross product
Vektor operator^(Vektor a, Vektor b) {
  return Vektor(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}
double udaljenost(Vektor a, Vektor b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  double dz = a.z - b.z;
  return sqrt(dx * dx + dy * dy + dz * dz);
}

double dot(Vektor a, Vektor b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double magnitude(Vektor v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }

void normaliziraj(Vektor &v) {
  double m = magnitude(v);
  if (m > 0) {
    v.x /= m;
    v.y /= m;
    v.z /= m;
  }
}
Vektor normalized(Vektor v) {
  Vektor rez;
  double m = magnitude(v);
  if (m > 0) {
    rez.x = v.x / m;
    rez.y = v.y / m;
    rez.z = v.z / m;
  }
  return rez;
}
Vektor closestLinePoint(Vektor point, Vektor a, Vektor b) {
  Vektor d = b - a;
  double t = ((point - a) * d) / (d * d); // projection along segment
  if (t < 0)
    t = 0;
  else if (t > 1)
    t = 1;
  return a + d * t;
}

double closestLineDistance(Vektor point, Vektor a, Vektor b) {
  Vektor closest = closestLinePoint(point, a, b);
  return magnitude(point - closest);
}
Vektor closestTrianglePoint(Vektor point, Vektor a, Vektor b, Vektor c) {

  Vektor n = normalized((b - a) ^ (c - a));
  double dist = (point - a) * n;
  Vektor proj = point - dist * n;

  Vektor v0 = b - a;
  Vektor v1 = c - a;
  Vektor v2 = proj - a;

  double d00 = v0 * v0;
  double d01 = v0 * v1;
  double d11 = v1 * v1;
  double d20 = v2 * v0;
  double d21 = v2 * v1;

  double denom = d00 * d11 - d01 * d01;
  double u = (d11 * d20 - d01 * d21) / denom;
  double v = (d00 * d21 - d01 * d20) / denom;
  double w = 1 - u - v;

  if (u >= 0 and v >= 0 and w >= 0 and u <= 1 and v <= 1 and w <= 1) {
    return proj;
  }

  Vektor pAB = closestLinePoint(point, a, b);
  Vektor pBC = closestLinePoint(point, b, c);
  Vektor pCA = closestLinePoint(point, c, a);

  double dAB = magnitude(point - pAB);
  double dBC = magnitude(point - pBC);
  double dCA = magnitude(point - pCA);

  if (dAB < dBC and dAB < dCA)
    return pAB;
  else if (dBC < dCA)
    return pBC;
  else
    return pCA;
}

class Kamera {
public:
  Vektor a;
  Vektor n;
  Vektor desno;
  Vektor gore;
  double fov;
  SDL_Renderer *renderer;

  void setN(Vektor v) {
    this->n.x = v.x;
    this->n.y = v.y;
    this->n.z = v.z;
    this->desno.x =
        this->n.y / sqrt(this->n.x * this->n.x + this->n.y * this->n.y);
    this->desno.y =
        -this->n.x / sqrt(this->n.x * this->n.x + this->n.y * this->n.y);
    this->desno.z = 0;
    this->gore = this->n ^ this->desno;

    normaliziraj(n);
    normaliziraj(gore);
    normaliziraj(desno);
  }

  Kamera(Vektor a, Vektor n, double fov, SDL_Renderer *renderer) {
    this->a = a;
    this->n = n;
    this->fov = fov;
    this->renderer = renderer;
  }

  bool clipLineAgainstNearPlane(Vektor &p1, Vektor &p2,
                                double nearPlane = 0.1) {
    if (p1.z < nearPlane && p2.z < nearPlane) {
      return false;
    }

    if (p1.z >= nearPlane && p2.z >= nearPlane) {
      return true;
    }

    Vektor *behind = (p1.z < nearPlane) ? &p1 : &p2;
    Vektor *inFront = (p1.z >= nearPlane) ? &p1 : &p2;

    double t = (nearPlane - behind->z) / (inFront->z - behind->z);

    behind->x = inFront->x * t + behind->x * (1 - t);
    behind->y = inFront->y * t + behind->y * (1 - t);
    behind->z = nearPlane;

    return true;
  }

  Tocka Tocka3d_To_Tocka2d(const Vektor &vertex) {
    Vektor toPoint = a - vertex;

    Vektor cameraSpace;
    cameraSpace.x = dot(toPoint, this->desno);
    cameraSpace.y = dot(toPoint, this->gore);
    cameraSpace.z = dot(toPoint, this->n);

    return Tocka(cameraSpace.x, cameraSpace.y, cameraSpace.z);
  }

  void spojiTocke2dClipped(const Tocka &a_cam, const Tocka &b_cam, int r, int g,
                           int bcol) {
    Vektor p1 = Vektor(a_cam.x, a_cam.y, a_cam.z);
    Vektor p2 = Vektor(b_cam.x, b_cam.y, b_cam.z);

    if (!clipLineAgainstNearPlane(p1, p2)) {
      return;
    }

    double scale1 = fov / p1.z;
    double screenX1 = p1.x * scale1;
    double screenY1 = p1.y * scale1;

    double scale2 = fov / p2.z;
    double screenX2 = p2.x * scale2;
    double screenY2 = p2.y * scale2;

    SDL_SetRenderDrawColor(renderer, r, g, bcol, 255);
    SDL_RenderDrawLine(renderer, 960 + (int)round(screenX1),
                       540 - (int)round(screenY1), 960 + (int)round(screenX2),
                       540 - (int)round(screenY2));
  }

  void rotate(double fi, double theta) {
    if (fi >= 360) {
      fi = 0.0001;
    }
    fi = fi * M_PI / 180.0;
    theta = theta * M_PI / 180.0;
    Vektor nigga;
    nigga.x = sin(theta) * cos(fi);
    nigga.y = sin(theta) * sin(fi);
    nigga.z = cos(theta);
    setN(nigga);
  }

  void move(double naprijed, double desno1, double gore1) {
    double mag = sqrt(n.x * n.x + n.y * n.y);
    a.x += naprijed * n.x / mag;
    a.y += naprijed * n.y / mag;

    a.x += desno1 * this->desno.x;
    a.y += desno1 * this->desno.y;
    a.z += desno1 * this->desno.z;

    a.z += gore1;
  }
};

class Graf {
public:
  double t = 0;
  double mass;
  vector<Vektor> vrhovi;
  vector<vector<int>> bridovi;
  vector<vector<int>> faces;

  Vektor velocity;

  Graf(double mass, double v_x, double v_y, double v_z) {
    this->mass = mass;
    this->velocity.x = v_x;
    this->velocity.y = v_y;
    this->velocity.z = v_z;
  }
  Graf() {
    this->velocity.x = 0;
    this->velocity.y = 0;
    this->velocity.z = 0;
  }

  virtual void shift(double x, double y, double z) {
    for (int i = 0; i < vrhovi.size(); i++) {
      vrhovi[i].x += x;
      vrhovi[i].y += y;
      vrhovi[i].z += z;
    }
  }

  void renderGraf(Kamera &kamera) {
    for (int i = 0; i < bridovi.size(); i++) {
      Tocka a = kamera.Tocka3d_To_Tocka2d(vrhovi[bridovi[i][0]]);
      Tocka b = kamera.Tocka3d_To_Tocka2d(vrhovi[bridovi[i][1]]);
      kamera.spojiTocke2dClipped(a, b, 255, 255, 255);
    }
  }
};

// made using chat gpt
class Kugla : public Graf {
public:
  int rings;   // latitude divisions (from pole to pole)
  int sectors; // longitude divisions (around)
  double radius;
  Vektor center;
  Vektor omega;
  double w;
  // vx,vy,vz - initial velocity (same style as Kocka)
  // radius - sphere radius
  // cx,cy,cz - center coordinates
  // rings - latitude subdivisions (>=2), sectors - longitude subdivisions (>=3)
  Kugla(Vektor omega, double w, double mass, double vx, double vy, double vz,
        double radius, double cx, double cy, double cz, int rings = 16,
        int sectors = 24) {
    this->omega = omega;
    this->w = w;
    this->mass = mass;
    this->velocity.x = vx;
    this->velocity.y = vy;
    this->velocity.z = vz;

    this->radius = radius;
    this->center.x = cx;
    this->center.y = cy;
    this->center.z = cz;
    this->rings = max(2, rings);
    this->sectors = max(3, sectors);

    // Create vertices: (rings + 1) latitude rows, each with sectors vertices
    // phi goes from 0..PI (north pole to south pole)
    for (int i = 0; i <= this->rings; ++i) {
      double phi = M_PI * (double)i / (double)this->rings; // 0..PI
      double sinPhi = sin(phi);
      double cosPhi = cos(phi);
      for (int j = 0; j < this->sectors; ++j) {
        double theta = 2.0 * M_PI * (double)j / (double)this->sectors; // 0..2PI
        double sinTheta = sin(theta);
        double cosTheta = cos(theta);

        double x = cx + radius * sinPhi * cosTheta;
        double y = cy + radius * sinPhi * sinTheta;
        double z = cz + radius * cosPhi;

        this->vrhovi.push_back(Vektor(x, y, z));
      }
    }

    // Create edges (bridovi)
    // For each vertex at (i,j):
    //  - connect to next longitude (i, j+1 mod sectors)
    //  - connect to next latitude (i+1, j) if not the last latitude row
    for (int i = 0; i <= this->rings; ++i) {
      for (int j = 0; j < this->sectors; ++j) {
        int idx = i * this->sectors + j;
        int nextLon = i * this->sectors + ((j + 1) % this->sectors);

        // longitude edge (wraps around)
        this->bridovi.push_back({idx, nextLon});

        // latitude edge (to the row below), if any
        if (i < this->rings) {
          int nextLat = (i + 1) * this->sectors + j;
          this->bridovi.push_back({idx, nextLat});
        }
      }
    }
  }
  void shift(double x, double y, double z) override {
    this->center.x += x;
    this->center.y += y;
    this->center.z += z;
    for (int i = 0; i < vrhovi.size(); i++) {
      vrhovi[i].x += x;
      vrhovi[i].y += y;
      vrhovi[i].z += z;
    }
  }
  void rotate_around_v() {
    normaliziraj(omega);

    for (int i = 0; i < this->vrhovi.size(); i++) {
      Vektor point = this->vrhovi[i];
      Vektor p_rel = point - center;
      this->vrhovi[i] = center + p_rel * cos(w) + (omega ^ p_rel) * sin(w) +
                        omega * (omega * p_rel) * (1 - cos(w));
    }
  }
};

class Scena {
public:
  vector<Graf> static_grafovi;
  vector<Kugla> moving_grafovi;

  int crosshair = 5;
  double g;
  double t = 0;
  Scena(double g) { this->g = g; }

  void addToStaticScene(Graf graf) { static_grafovi.push_back(graf); }
  void addToMovingScene(Kugla graf) { moving_grafovi.push_back(graf); }

  void renderScene(Kamera &kamera, int delay) {
    for (int i = 0; i < static_grafovi.size(); i++) {
      static_grafovi[i].renderGraf(kamera);
    }

    for (int i = 0; i < moving_grafovi.size(); i++) {
      moving_grafovi[i].renderGraf(kamera);
    }

    SDL_RenderDrawLine(kamera.renderer, 960 - crosshair, 540, 960 + crosshair,
                       540);
    SDL_RenderDrawLine(kamera.renderer, 960, 540 - crosshair, 960,
                       540 + crosshair);

    SDL_RenderPresent(kamera.renderer);
    SDL_Delay(delay);
    SDL_SetRenderDrawColor(kamera.renderer, 0, 0, 0, 255);
    SDL_RenderClear(kamera.renderer);
  }
  void nacrtajKordinate(Kamera &kamera) {
    Vektor ishodiste = Vektor(0, 0, 0);
    Vektor x = Vektor(5000, 0, 0);
    Vektor y = Vektor(0, 5000, 0);
    Vektor z = Vektor(0, 0, 5000);

    Tocka a = kamera.Tocka3d_To_Tocka2d(ishodiste);
    Tocka b = kamera.Tocka3d_To_Tocka2d(x);
    kamera.spojiTocke2dClipped(a, b, 255, 0, 0);

    b = kamera.Tocka3d_To_Tocka2d(y);
    kamera.spojiTocke2dClipped(a, b, 0, 255, 0);

    b = kamera.Tocka3d_To_Tocka2d(z);
    kamera.spojiTocke2dClipped(a, b, 0, 0, 255);
  }

  void animateScene(double deltaTime) {
    const double restitution = 0.6;

    // movement
    for (int i = 0; i < moving_grafovi.size(); i++) {
      moving_grafovi[i].rotate_around_v();
      moving_grafovi[i].shift(moving_grafovi[i].velocity.x * deltaTime * 60,
                              moving_grafovi[i].velocity.y * deltaTime * 60,
                              moving_grafovi[i].velocity.z * deltaTime * 60);

      // gravity
      moving_grafovi[i].velocity.z -= g * deltaTime;
    }

    for (int iter = 0; iter < 5; iter++) { // multiple iterations for stability
      for (int i = 0; i < moving_grafovi.size(); i++) {
        for (int j = i + 1; j < moving_grafovi.size(); j++) {
          Vektor diff = moving_grafovi[i].center - moving_grafovi[j].center;
          double dist = magnitude(diff);
          double rsum = moving_grafovi[i].radius + moving_grafovi[j].radius;

          if (dist < rsum && dist > 0) {
            Vektor n = (1 / dist) * diff;

            Vektor v1 = moving_grafovi[i].velocity;
            Vektor v2 = moving_grafovi[j].velocity;
            double m1 = moving_grafovi[i].mass;
            double m2 = moving_grafovi[j].mass;

            double vrel = (v1 - v2) * n;
            if (vrel < 0) {
              double jImpulse =
                  -(1.0 + restitution) * vrel / (1.0 / m1 + 1.0 / m2);
              moving_grafovi[i].velocity = v1 + (jImpulse / m1) * n;
              moving_grafovi[j].velocity = v2 - (jImpulse / m2) * n;
            }

            // position correction with slop
            const double percent = 0.8; // solve 80% this frame
            const double slop = 0.01;
            double penetration = rsum - dist;
            if (penetration > slop) {
              penetration = (penetration - slop) * percent;
              Vektor correction = (penetration / (m1 + m2)) * n;
              moving_grafovi[i].shift(correction.x * m2, correction.y * m2,
                                      correction.z * m2);
              moving_grafovi[j].shift(-correction.x * m1, -correction.y * m1,
                                      -correction.z * m1);
            }
          }
        }
      }
    }

    // collision with static objects
    for (int i = 0; i < moving_grafovi.size(); i++) {

      for (int j = 0; j < static_grafovi.size(); j++) {
        for (int k = 0; k < static_grafovi[j].faces.size(); k++) {
          vector<int> trojka = static_grafovi[j].faces[k];

          Vektor t1 = static_grafovi[j].vrhovi[trojka[0]];
          Vektor t2 = static_grafovi[j].vrhovi[trojka[1]];
          Vektor t3 = static_grafovi[j].vrhovi[trojka[2]];
          Vektor impactPoint =
              closestTrianglePoint(moving_grafovi[i].center, t1, t2, t3);
          double pen = udaljenost(impactPoint, moving_grafovi[i].center) -
                       moving_grafovi[i].radius;
          if (pen < 0) {

            Vektor n = normalized(impactPoint - moving_grafovi[i].center);
            Vektor vel = moving_grafovi[i].velocity;
            Vektor korekcija1 = pen * n;
            moving_grafovi[i].shift(korekcija1.x, korekcija1.y, korekcija1.z);
            moving_grafovi[i].velocity =
                vel - (1 + restitution) * (vel * n) * n;
          }
        }
      }
    }
  }
};

class Kocka : public Graf {
public:
  int velicina;
  double xsred;
  double ysred;
  double zsred;

  Kocka(double v_x, double v_y, double v_z, double velicina, double xsred,
        double ysred, double zsred) {
    this->velocity.x = v_x;
    this->velocity.y = v_y;
    this->velocity.z = v_z;
    int a[8][3] = {{1, 1, 1},  {1, 1, -1},  {1, -1, 1},  {1, -1, -1},
                   {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};
    for (int i = 0; i < 8; i++) {
      this->vrhovi.push_back(Vektor(xsred + a[i][0] * velicina / 2,
                                    ysred + a[i][1] * velicina / 2,
                                    zsred + a[i][2] * velicina / 2));
    }
    this->bridovi = {{0, 1}, {1, 3}, {3, 2}, {2, 0}, {4, 5}, {5, 7},
                     {7, 6}, {6, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

    this->faces = {/*{0, 2, 6}, {0, 6, 4},*/ {1, 5, 7},
                   {1, 7, 3},
                   {4, 6, 7},
                   {4, 7, 5},
                   {0, 1, 3},
                   {0, 3, 2},
                   {0, 4, 5},
                   {0, 5, 1},
                   {2, 3, 7},
                   {2, 7, 6}};
  }
};

class Kocka1 : public Graf {
public:
  int velicina;
  double xsred;
  double ysred;
  double zsred;

  Kocka1(double v_x, double v_y, double v_z, double velicina, double xsred,
         double ysred, double zsred) {
    this->velocity.x = v_x;
    this->velocity.y = v_y;
    this->velocity.z = v_z;
    int a[8][3] = {{1, 1, 3},  {1, 1, -1},  {1, -1, 3},  {1, -1, -1},
                   {-1, 1, 3}, {-1, 1, -1}, {-1, -1, 3}, {-1, -1, -1}};
    for (int i = 0; i < 8; i++) {
      this->vrhovi.push_back(Vektor(xsred + a[i][0] * velicina / 2,
                                    ysred + a[i][1] * velicina / 2,
                                    zsred + a[i][2] * velicina / 2));
    }
    this->bridovi = {{0, 1}, {1, 3}, {3, 2}, {2, 0}, {4, 5}, {5, 7},
                     {7, 6}, {6, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

    this->faces = {{1, 5, 7}, {1, 7, 3}, {4, 6, 7}, {4, 7, 5}, {0, 1, 3},
                   {0, 3, 2}, {0, 4, 5}, {0, 5, 1}, {2, 3, 7}, {2, 7, 6}};
  }
};

int main() {
  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window *window = SDL_CreateWindow("Game", SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED, 1920, 1080, 0);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);

  SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
  SDL_RenderClear(renderer);

  Kamera kamera =
      Kamera(Vektor(684.724, -2781.49, 1880), Vektor(-1, 0, 0), 900, renderer);

  bool running = true;
  bool freeze = false;
  int sprinting = 1;
  SDL_Event event;
  double rotationSpeed = 1;
  Scena scena(7 * 9.81);
  int raspon = 100000;

  int previosMouseY = 1080 / 2;
  int previosMouseX = 1920 / 2;
  double fi = -82.6;
  double theta = 62.8;
  double brzina = 1;
  double brzina1 = 20;
  SDL_SetRelativeMouseMode(SDL_TRUE);
  int deltaX;
  int deltaY;

  double t = 0;
  double velocity = 0;
  double pod = 0;
  int f = 0;
  Kocka kocka1(0, 0, 0, 4000, 0, 0, 0);
  scena.addToStaticScene(kocka1);

  /*Kugla kugla1(Vektor(0, 0, 1), 0, 10000, 0, 0, 0, 800, 0, 0, 900, 7, 13);
  scena.addToMovingScene(kugla1);*/
  /* Povrsina p1 = Povrsina(Vektor(0.5, 0, 1), Vektor(0, 0, -5000));
   Povrsina p2 = Povrsina(Vektor(0, 0, -1), Vektor(0, 0, 5000));
   Povrsina p3 = Povrsina(Vektor(0, 1, 0), Vektor(0, -5000, 0));
   Povrsina p4 = Povrsina(Vektor(0, -1, 0), Vektor(0, 5000, 0));
   Povrsina p5 = Povrsina(Vektor(1, 0, 0), Vektor(-5000, 0, 0));
   Povrsina p6 = Povrsina(Vektor(-1, 0, 0), Vektor(5000, 0, 0));

   scena.addToPovrsine(p1);
   scena.addToPovrsine(p2);
   scena.addToPovrsine(p3);
   scena.addToPovrsine(p4);
   scena.addToPovrsine(p5);
   scena.addToPovrsine(p6);*/

  Uint64 NOW = SDL_GetPerformanceCounter();
  Uint64 LAST = 0;
  double deltaTime = 0;

  while (running) {
    /*if (f == 4) {
      Kugla kugla(Vektor(0, 0, 1), 0, 10, 10, 0, 0, 200, 0, 0, 4000, 5, 10);
      scena.addToMovingScene(kugla);
      f = 0;
    }
    f++;*/
    while (SDL_PollEvent(&event)) {

      if (event.type == SDL_QUIT) {
        running = false;
      } else if (event.type == SDL_KEYDOWN &&
                 event.key.keysym.sym == SDLK_ESCAPE) {

        running = false;

      } else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_p) {
        if (freeze == true)
          freeze = false;
        else {
          freeze = true;
        }
      } else if (event.type == SDL_MOUSEBUTTONDOWN &&
                 event.button.button == SDL_BUTTON_LEFT) {
        Kugla kugla(Vektor(0, 0, 1), 0, 100, -30 * kamera.n.x, -30 * kamera.n.y,
                    -30 * kamera.n.z, 200, kamera.a.x, kamera.a.y, kamera.a.z,
                    5, 10);
        scena.addToMovingScene(kugla);
      } else if (event.type == SDL_MOUSEBUTTONDOWN &&
                 event.button.button == SDL_BUTTON_RIGHT) {
        Kugla kugla(Vektor(0, 0, 1), 0, 100, -120 * kamera.n.x,
                    -120 * kamera.n.y, -120 * kamera.n.z, 200, kamera.a.x,
                    kamera.a.y, kamera.a.z, 5, 10);
        scena.addToMovingScene(kugla);
      }
    }

    const Uint8 *state = SDL_GetKeyboardState(NULL);

    if (state[SDL_SCANCODE_LCTRL]) {
      sprinting = 5;
    } else {
      sprinting = 1;
    }
    if (state[SDL_SCANCODE_W])
      kamera.move(-brzina1 * sprinting, 0, 0);
    if (state[SDL_SCANCODE_S])
      kamera.move(brzina1 * sprinting, 0, 0);
    if (state[SDL_SCANCODE_A])
      kamera.move(0, brzina1 * sprinting, 0);
    if (state[SDL_SCANCODE_D])
      kamera.move(0, -brzina1 * sprinting, 0);
    if (state[SDL_SCANCODE_LSHIFT])
      kamera.move(0, 0, -brzina1 * sprinting);
    if (state[SDL_SCANCODE_SPACE])
      kamera.move(0, 0, brzina1 * sprinting);

    if (event.type == SDL_MOUSEMOTION) {
      SDL_GetRelativeMouseState(&deltaX, &deltaY);
      theta -= (double)deltaY / 10;
      fi -= (double)deltaX / 10;

      if (theta > 179.9999) {
        theta = 179.9999;
      } else if (theta < 0.0001) {
        theta = 0.0001;
      }
    }

    fi = fmod(fi, 360);
    if (fi == 0 or fi == 360) {
      fi = 0.0001;
    }
    if (theta == 0 or theta == 180) {
      theta = 0.0001;
    }
    kamera.rotate(fi, theta);
    // cout << fi << " " << theta << endl;
    // kamera.a.printVrh();

    scena.nacrtajKordinate(kamera);
    LAST = NOW;
    NOW = SDL_GetPerformanceCounter();
    deltaTime =
        (double)((NOW - LAST) * 1000 / (double)SDL_GetPerformanceFrequency());
    deltaTime /= 1000;

    if (freeze == false) {
      scena.animateScene(deltaTime);
    }
    scena.renderScene(kamera, 16);
  }

  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
