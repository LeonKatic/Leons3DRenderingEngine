#include <SDL2/SDL.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_scancode.h>
#include <SDL2/SDL_stdinc.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <math.h>

#include <stdlib.h>
#include <unistd.h>
#include <vector>
using namespace std;

int idGlobal = 0;
class Vektor;
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

class Vrh {
public:
  double x;
  double y;
  double z;
  Vrh() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
  }
  Vrh(double xkor, double ykor, double zkor) {
    this->x = xkor;
    this->y = ykor;
    this->z = zkor;
  }

  Vrh(const class Vektor &vektor);
  void printVrh() {
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
  Vektor(const Vrh &vrh) {
    x = vrh.x;
    y = vrh.y;
    z = vrh.z;
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
Vrh::Vrh(const Vektor &vektor) {
  this->x = vektor.x;
  this->y = vektor.y;
  this->z = vektor.z;
}
Vektor oduzmi(Vrh a, Vrh b) {
  Vektor rez = Vektor();
  rez.x = b.x - a.x;
  rez.y = b.y - a.y;
  rez.z = b.z - a.z;
  return rez;
}
Vektor oduzmiVek(Vektor a, Vektor b) {
  Vektor rez = Vektor();
  rez.x = b.x - a.x;
  rez.y = b.y - a.y;
  rez.z = b.z - a.z;
  return rez;
}
Vrh zbroji(Vrh a, Vrh b) {
  Vrh rez = Vrh();
  rez.x = b.x + a.x;
  rez.y = b.y + a.y;
  rez.z = b.z + a.z;
  return rez;
}
/*Vektor zbrojivek(Vektor a, Vektor b) {
  Vektor rez = Vektor();
  rez.x = b.x + a.x;
  rez.y = b.y + a.y;
  rez.z = b.z + a.z;
  return rez;
}*/
Vrh skalar(double a, Vrh vrh) {
  Vrh rez = Vrh();
  rez.x = a * vrh.x;
  rez.y = a * vrh.y;
  rez.z = a * vrh.z;
  return rez;
}
Vektor skalarVek(double a, Vektor v) {
  Vektor rez = Vektor();
  rez.x = a * v.x;
  rez.y = a * v.y;
  rez.z = a * v.z;
  return rez;
}

Vektor cross(Vektor a, Vektor b) {
  Vektor rez;
  rez.x = a.y * b.z - a.z * b.y;
  rez.y = a.z * b.x - a.x * b.z;
  rez.z = a.x * b.y - a.y * b.x;
  return rez;
}
double udaljenost(Vrh a, Vrh b) {
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

class Kamera {
public:
  Vrh a;
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
    this->gore = cross(this->n, this->desno);

    normaliziraj(n);
    normaliziraj(gore);
    normaliziraj(desno);
  }

  Kamera(Vrh a, Vektor n, double fov, SDL_Renderer *renderer) {
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

  Tocka Tocka3d_To_Tocka2d(const Vrh &vertex) {
    Vektor toPoint = oduzmi(vertex, a);

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
  vector<Vrh> vrhovi;
  vector<pair<int, int>> bridovi;
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

  Vrh minVrh() {
    Vrh min = this->vrhovi[0];
    for (int i = 0; i < this->vrhovi.size(); i++) {
      if (this->vrhovi[i].z < min.z)
        min = this->vrhovi[i];
    }
    return min;
  }

  void renderGraf(Kamera &kamera) {
    for (int i = 0; i < bridovi.size(); i++) {
      Tocka a = kamera.Tocka3d_To_Tocka2d(vrhovi[bridovi[i].first]);
      Tocka b = kamera.Tocka3d_To_Tocka2d(vrhovi[bridovi[i].second]);
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
  Vrh center;

  // vx,vy,vz - initial velocity (same style as Kocka)
  // radius - sphere radius
  // cx,cy,cz - center coordinates
  // rings - latitude subdivisions (>=2), sectors - longitude subdivisions (>=3)
  Kugla(double mass, double vx, double vy, double vz, double radius, double cx,
        double cy, double cz, int rings = 16, int sectors = 24) {
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

        this->vrhovi.push_back(Vrh(x, y, z));
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
};

class Scena {
public:
  vector<Graf> static_grafovi;
  vector<Kugla> moving_grafovi;

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
    SDL_RenderPresent(kamera.renderer);
    SDL_Delay(delay);
    SDL_SetRenderDrawColor(kamera.renderer, 0, 0, 0, 255);
    SDL_RenderClear(kamera.renderer);
  }
  void nacrtajKordinate(Kamera &kamera) {
    Vrh ishodiste = Vrh(0, 0, 0);
    Vrh x = Vrh(5000, 0, 0);
    Vrh y = Vrh(0, 5000, 0);
    Vrh z = Vrh(0, 0, 5000);

    Tocka a = kamera.Tocka3d_To_Tocka2d(ishodiste);
    Tocka b = kamera.Tocka3d_To_Tocka2d(x);
    kamera.spojiTocke2dClipped(a, b, 255, 0, 0);

    b = kamera.Tocka3d_To_Tocka2d(y);
    kamera.spojiTocke2dClipped(a, b, 0, 255, 0);

    b = kamera.Tocka3d_To_Tocka2d(z);
    kamera.spojiTocke2dClipped(a, b, 0, 0, 255);
  }
  void animateScene(double deltaTime) {

    for (int i = 0; i < moving_grafovi.size(); i++) {

      moving_grafovi[i].velocity.printVektor();
      // moving_grafovi[i].center.printVrh();
      moving_grafovi[i].shift(moving_grafovi[i].velocity.x * deltaTime * 60,
                              moving_grafovi[i].velocity.y * deltaTime * 60,
                              moving_grafovi[i].velocity.z * deltaTime * 60);
      // force of gravity
      moving_grafovi[i].velocity.z -= g * deltaTime;
    }
    // detect collision

    for (int i = 0; i < moving_grafovi.size(); i++) {
      for (int j = i + 1; j < moving_grafovi.size(); j++) {
        double dist =
            udaljenost(moving_grafovi[i].center, moving_grafovi[j].center);
        double delta =
            dist - (moving_grafovi[i].radius + moving_grafovi[j].radius);
        if (delta < -1e-6) {
          double m1 = moving_grafovi[i].mass;
          double m2 = moving_grafovi[j].mass;
          double m_uk = m1 + m2;

          Vektor n_col =
              oduzmi(moving_grafovi[j].center, moving_grafovi[i].center);
          normaliziraj(n_col);

          Vektor v1 = moving_grafovi[i].velocity;
          Vektor v2 = moving_grafovi[j].velocity;

          // elastic collision
          double v1n = (v1 * n_col);
          double v2n = (v2 * n_col);

          double p = 2.0 * (v1n - v2n) / (m1 + m2);

          moving_grafovi[i].velocity = v1 - n_col * (p * m2);
          moving_grafovi[j].velocity = v2 + n_col * (p * m1);

          double correctionFactor = 0.5;
          Vektor correction = n_col * (abs(delta) * correctionFactor);

          moving_grafovi[i].shift(-correction.x * (m2 / m_uk),
                                  -correction.y * (m2 / m_uk),
                                  -correction.z * (m2 / m_uk));
          moving_grafovi[j].shift(correction.x * (m1 / m_uk),
                                  correction.y * (m1 / m_uk),
                                  correction.z * (m1 / m_uk));
        }
      }
    }

    // collsion on plane z=0;
    for (int i = 0; i < moving_grafovi.size(); i++) {
      if (moving_grafovi[i].minVrh().z <= 0) {
        moving_grafovi[i].velocity.z *= -0.9;
      }
    }
  }
  void animateSceneCHATGPT(double deltaTime) {
    const double restitution =
        0.9;                    // 1.0 = perfectly elastic, <1 = damped bounce
    const double percent = 0.2; // positional correction strength (20%)
    const double slop = 0.01;   // penetration tolerance

    // --- integrate motion (gravity + movement) ---
    for (int i = 0; i < moving_grafovi.size(); i++) {
      moving_grafovi[i].shift(moving_grafovi[i].velocity.x * deltaTime * 60,
                              moving_grafovi[i].velocity.y * deltaTime * 60,
                              moving_grafovi[i].velocity.z * deltaTime * 60);

      // gravity
      moving_grafovi[i].velocity.z -= g * deltaTime;
    }

    // --- sphere-sphere collisions ---
    for (int i = 0; i < moving_grafovi.size(); i++) {
      for (int j = i + 1; j < moving_grafovi.size(); j++) {
        Vektor n = oduzmi(moving_grafovi[j].center, moving_grafovi[i].center);
        double dist = magnitude(n);
        double rsum = moving_grafovi[i].radius + moving_grafovi[j].radius;

        if (dist < rsum && dist > 1e-8) {
          normaliziraj(n);

          double m1 = moving_grafovi[i].mass;
          double m2 = moving_grafovi[j].mass;
          double invM1 = 1.0 / m1;
          double invM2 = 1.0 / m2;

          // relative velocity
          Vektor rv = moving_grafovi[i].velocity - moving_grafovi[j].velocity;
          double velAlongNormal = rv * n;

          // only resolve if moving toward each other
          if (velAlongNormal < 0) {
            double j_imp = -(1 + restitution) * velAlongNormal;
            j_imp /= (invM1 + invM2);

            Vektor impulse = n * j_imp;
            moving_grafovi[i].velocity =
                moving_grafovi[i].velocity + impulse * invM1;
            moving_grafovi[j].velocity =
                moving_grafovi[j].velocity - impulse * invM2;
          }

          // positional correction (Baumgarte stabilization)
          double penetration = rsum - dist;
          double correctionMag =
              max(penetration - slop, 0.0) / (invM1 + invM2);
          Vektor correction = n * (correctionMag * percent);

          moving_grafovi[i].shift(-correction.x * invM1, -correction.y * invM1,
                                  -correction.z * invM1);
          moving_grafovi[j].shift(correction.x * invM2, correction.y * invM2,
                                  correction.z * invM2);
        }
      }
    }

    // --- collisions with floor (z = 0) ---
    for (int i = 0; i < moving_grafovi.size(); i++) {
      if (moving_grafovi[i].minVrh().z <= 0) {
        double pen = 0.0 - moving_grafovi[i].minVrh().z;
        if (pen > 0.0) {
          moving_grafovi[i].shift(0, 0, pen + 0.001);
        }
        moving_grafovi[i].velocity.z *= -restitution;
      }
    }

    // --- debug log ---
    for (int i = 0; i < moving_grafovi.size(); i++) {
      cout << i << " ";
      moving_grafovi[i].velocity.printVektor();
      moving_grafovi[i].t += 0.01;
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
      this->vrhovi.push_back(Vrh(xsred + a[i][0] * velicina / 2,
                                 ysred + a[i][1] * velicina / 2,
                                 zsred + a[i][2] * velicina / 2));
    }
    this->bridovi = {{0, 1}, {1, 3}, {3, 2}, {2, 0}, {4, 5}, {5, 7},
                     {7, 6}, {6, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
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
      Kamera(Vrh(684.724, -2781.49, 1880), Vektor(-1, 0, 0), 900, renderer);

  bool running = true;
  bool freeze = false;
  int sprinting = 1;
  SDL_Event event;
  double rotationSpeed = 1;
  Scena scena(9.81);
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
  Kocka kocka1(0, 0, 0, 10000, 0, 0, 0);
  scena.addToStaticScene(kocka1);
  Kugla kugla1(10000, 0, 0, 0, 800, 0, 0, 900, 7, 13);
  scena.addToMovingScene(kugla1);
  /*
    Kugla kugla2(200, 30, 0, 0, 500, -2100, -40, 500, 5, 10);
    scena.addToMovingScene(kugla2);

    Kugla kugla3(200, 30, 0, 0, 500, 0, 0, 0, 5, 10);
    scena.addToMovingScene(kugla3);*/
  Uint64 NOW = SDL_GetPerformanceCounter();
  Uint64 LAST = 0;
  double deltaTime = 0;
  while (running) {
    /*if (velocity < 1000) {
      velocity += (double)1 / 5;
    }

    if (scena.grafovi[0].minVrh().z < pod) {
      scena.grafovi[0].shift(0, 0, velocity);
    }
    if (scena.grafovi[0].minVrh().z > pod) {
      scena.grafovi[0].shift(0, 0, pod - scena.grafovi[0].minVrh().z);
    }

    scena.grafovi[0].shift(5, 0, 0);
    t++;
    if (freeze == false) {
      f++;
      if (f > 50) {
        cout << f << endl;
        Kugla kugla1(100, -30, 0, 0, 200, 2100, 0, 500, 5, 10);
        scena.addToMovingScene(kugla1);

        Kugla kugla1(rand()%400,rand() % 90 - 45, rand() % 90 - 45, rand() %
    100, rand() % 200+100, rand() % 900 - 450, rand() % 900 - 450, rand() % 900
    - 450, 5, 10); // vx,vy,vz, radius, cx,cy,cz, rings, sectors
        scena.addToMovingScene(kugla1);

        f = 0;
      }
    }*/

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
        Kugla kugla(1000, -30 * kamera.n.x, -30 * kamera.n.y, -30 * kamera.n.z,
                    200, kamera.a.x, kamera.a.y, kamera.a.z, 5, 10);
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
      scena.animateSceneCHATGPT(deltaTime);
    }
    scena.renderScene(kamera, 16);
  }

  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}