#include <array>
#include <cmath>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <cstdio>  //        Remove "-fopenmp" for g++ version < 4.2
#include <cstdlib> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <fstream>
#include <iostream>
#include <random>
#include <span>
#include <vector>
struct __attribute__((aligned(32))) Vec
{                   // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    explicit Vec(double x_ = 0, double y_ = 0, double z_ = 0) noexcept : x(x_), y(y_), z(z_){};
    [[nodiscard]] Vec operator+(const Vec& b) const noexcept { return Vec(x + b.x, y + b.y, z + b.z); }
    [[nodiscard]] Vec operator-(const Vec& b) const noexcept { return Vec(x - b.x, y - b.y, z - b.z); }
    [[nodiscard]] Vec operator*(double b) const noexcept { return Vec(x * b, y * b, z * b); }
    [[nodiscard]] Vec mult(const Vec& b) const noexcept { return Vec(x * b.x, y * b.y, z * b.z); }
    [[nodiscard]] Vec& norm() noexcept { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    [[nodiscard]] double dot(const Vec& b) const noexcept { return x * b.x + y * b.y + z * b.z; } // cross:
    [[nodiscard]] Vec operator%(Vec& b) const noexcept { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct __attribute__((aligned(64))) Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()
struct __attribute__((aligned(128))) Sphere
{
    Vec p, e, c; // position, emission, color
    double rad;  // radius
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(Vec p_, Vec e_, Vec c_, double rad_, Refl_t refl_) noexcept : p(p_), e(e_), c(c_), rad(rad_), refl(refl_) {}
    [[nodiscard]] double intersect(const Ray& r) const
    {                     // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double eps = 1e-4;
        double b = op.dot(r.d);
        double det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
        {
            return 0;
        }
        det = sqrt(det);
        double t = b - det;
        return t > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};
const std::array<Sphere, 9> spheres = {
    // Scene: radius, position, emission, color, material
    Sphere(Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), 1e5, DIFF),   // Left
    Sphere(Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), 1e5, DIFF), // Rght
    Sphere(Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), 1e5, DIFF),         // Back
    Sphere(Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), 1e5, DIFF),               // Frnt
    Sphere(Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), 1e5, DIFF),         // Botm
    Sphere(Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), 1e5, DIFF), // Top
    Sphere(Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, 16.5, SPEC),        // Mirr
    Sphere(Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, 16.5, REFR),        // Glas
    Sphere(Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), 600, DIFF)     // Lite
};
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }
inline int toInt(double x) { return static_cast<int>(lround(pow(clamp(x), 1 / 2.2) * 255)); }
inline bool intersect(const Ray& r, double& t, int& id)
{
    double n = static_cast<double>(sizeof(spheres)) / sizeof(Sphere);
    double inf = t = 1e20;
    for (int i = static_cast<int>(n); i >= 0; i--)
    {
        double d = spheres[static_cast<uint32_t>(i)].intersect(r);
        if (d > 0.0 && d < t)
        {
            t = d;
            id = i;
        }
    }
    return t < inf;
}
Vec radiance(const Ray& r, int depth, std::mt19937& gen)
{
    std::uniform_real_distribution<> dis(0.0F, 1.0F);
    double t = 1e20; // distance to intersection
    int id = 0;      // id of intersected object
    if (!intersect(r, t, id))
    {
        return Vec(); // if miss, return black
    }
    const Sphere& obj = spheres[static_cast<uint32_t>(id)]; // the hit object
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    Vec f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                        : f.z; // max refl
    if (++depth > 5)
    {
        if (dis(gen) < p)
        {
            f = f * (1 / p);
        } else
        {
            return obj.e; // R.R.
        }
    }
    if (obj.refl == DIFF)
    { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * dis(gen);
        double r2 = dis(gen);
        double r2s = sqrt(r2);
        Vec w = nl;
        Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
        Vec v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, gen));
    }
    if (obj.refl == SPEC) // Ideal SPECULAR reflection
    {
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, gen));
    }
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1;
    double nt = 1.5;
    double nnt = into ? nc / nt : nt / nc;
    double ddn = r.d.dot(nl);
    double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    if (cos2t < 0) // Total internal reflection
    {
        return obj.e + f.mult(radiance(reflRay, depth, gen));
    }
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc;
    double b = nt + nc;
    double R0 = a * a / (b * b);
    double c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c;
    double Tr = 1 - Re;
    double P = .25 + .5 * Re;
    double RP = Re / P;
    double TP = Tr / (1 - P);
    return obj.e + f.mult(depth > 2 ? (dis(gen) < P ? // Russian roulette
                                           radiance(reflRay, depth, gen) * RP
                                                    : radiance(Ray(x, tdir), depth, gen) * TP)
                                    : radiance(reflRay, depth, gen) * Re + radiance(Ray(x, tdir), depth, gen) * Tr);
}

inline void ProcessPixel(uint32_t h, uint32_t y, uint32_t w,
                         uint32_t x, int samps, std::uniform_real_distribution<>& dis,
                         const Vec& cx, const Vec& cy, const Ray& cam, std::mt19937& gen,
                         std::vector<Vec>& c)
{
    for (uint32_t sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
    {
        for (int sx = 0; sx < 2; sx++)
        {
            // 2x2 subpixel cols
            Vec r;
            for (int s = 0; s < samps; s++)
            {
                double r1 = 2 * dis(gen);
                double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                double r2 = 2 * dis(gen);
                double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, gen) * (1. / samps);
            } // Camera rays are pushed ^^^^^ forward to start in interior
            c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
        }
    }
}

int main(int argc, char* argv[])
{
    constexpr uint32_t w = 1024;
    constexpr uint32_t h = 768;
    int samps = argc == 2 ? std::stoi(std::span(argv, static_cast<uint32_t>(argc))[1]) / 4 : 1; // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());                                  // cam pos, dir
    Vec cx = Vec(w * .5135 / h);
    Vec cy = (cx % cam.d).norm() * .5135;
    // Vec r;
    std::vector<Vec> c(static_cast<size_t>(w * h)); // = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1)       // OpenMP
    for (uint32_t y = 0; y < h; y++)
    { // Loop over image rows
        std::cerr << "\rRendering (" << samps * 4 << " spp) " << 100. * y / (h - 1);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0F, 1.0F);
        for (uint32_t x = 0; x < w; x++) // Loop cols
        {
            ProcessPixel(h, y, w, x, samps, dis, cx, cy, cam, gen, c);
        }
    }
    std::ofstream f;
    f.open("image.ppm"); // Write image to PPM file.
    f << "P3\n"
      << w << " " << h << "\n"
      << 255 << "\n";
    for (auto printVec : c)
    {
        f << toInt(printVec.x) << " " << toInt(printVec.y) << " " << toInt(printVec.z) << " ";
    }
}
