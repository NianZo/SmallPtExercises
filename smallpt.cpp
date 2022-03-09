#include <array>
#include <cmath>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <cstdio>  //        Remove "-fopenmp" for g++ version < 4.2
#include <cstdlib> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <fstream>
#include <iostream>
#include <random>
#include <span>
#include <sstream>
#include <vector>

constexpr auto seed()
{
  std::uint64_t shifted = 0;

  for( const auto c : __TIME__ )
  {
    shifted <<= 8;
    shifted |= c;
  }

  return shifted;
}

struct PCG
{
  struct pcg32_random_t { std::uint64_t state=0;  std::uint64_t inc=seed(); };
  pcg32_random_t rng;
  typedef std::uint32_t result_type;

  constexpr result_type operator()()
  {
    return pcg32_random_r();
  }

  static result_type constexpr min()
  {
    return std::numeric_limits<result_type>::min();
  }

  static result_type constexpr max()
  {
    return std::numeric_limits<result_type>::min();
  }

  private:
  constexpr std::uint32_t pcg32_random_r()
  {
    std::uint64_t oldstate = rng.state;
    // Advance internal state
    rng.state = oldstate * 6364136223846793005ULL + (rng.inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    std::uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  }

};

constexpr auto get_random(int count)
{
  PCG pcg;
  while(count > 0){
    pcg();
    --count;
  }
  return static_cast<double>(pcg()) / std::numeric_limits<uint32_t>::max();
}

struct __attribute__((aligned(32))) Vec
{                   // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    explicit constexpr Vec(double x_ = 0, double y_ = 0, double z_ = 0) noexcept : x(x_), y(y_), z(z_){};
    [[nodiscard]] constexpr Vec operator+(const Vec& b) const noexcept { return Vec(x + b.x, y + b.y, z + b.z); }
    [[nodiscard]] constexpr Vec operator-(const Vec& b) const noexcept { return Vec(x - b.x, y - b.y, z - b.z); }
    [[nodiscard]] constexpr Vec operator*(double b) const noexcept { return Vec(x * b, y * b, z * b); }
    [[nodiscard]] constexpr Vec mult(const Vec& b) const noexcept { return Vec(x * b.x, y * b.y, z * b.z); }
    [[nodiscard]] constexpr Vec norm() const noexcept { return Vec(*this) * (1 / sqrt(x * x + y * y + z * z)); }
    [[nodiscard]] constexpr double dot(const Vec& b) const noexcept { return x * b.x + y * b.y + z * b.z; } // cross:
    [[nodiscard]] constexpr Vec operator%(const Vec& b) const noexcept { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct __attribute__((aligned(64))) Ray
{
    const Vec o, d;
    constexpr Ray(const Vec& o_, const Vec& d_) : o(o_), d(d_) {}
};
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()
struct __attribute__((aligned(128))) Sphere
{
    const Vec p, e, c; // position, emission, color
    const double rad;  // radius
    const Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    constexpr Sphere(const Vec& p_, const Vec& e_, const Vec& c_, double rad_, Refl_t refl_) noexcept : p(p_), e(e_), c(c_), rad(rad_), refl(refl_) {}
    [[nodiscard]] constexpr double intersect(const Ray& r) const
    {                           // returns distance, 0 if nohit
        const Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        constexpr double eps = 1e-4;
        const double b = op.dot(r.d);
        const double det2 = b * b - op.dot(op) + rad * rad;
        if (det2 < 0)
        {
            return 0;
        }
        const double det = sqrt(det2);
        const double t1 = b - det;
        const double t2 = b + det;
        return t1 > eps ? t1 : (t2 > eps ? t2 : 0);
    }
};
constexpr std::array<Sphere, 9> spheres = {
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
inline constexpr double clamp(const double x) { return x < 0 ? 0 : x > 1 ? 1
                                                               : x; }
inline constexpr int toInt(const double x) { return static_cast<int>(lround(pow(clamp(x), 1 / 2.2) * 255)); }
struct Intersection
{
	double t;
	uint32_t id;
	bool hit;
};
inline constexpr Intersection intersect(const Ray& r)
{
    //constexpr double n = static_cast<double>(sizeof(spheres)) / sizeof(Sphere);
    constexpr double inf = 1e20;
    double t = 1e20;
    uint32_t id = 0;
    for (uint32_t i = 0; i < spheres.size(); i++)
    {
        const double d = spheres[i].intersect(r);
        if (d > 0.0 && d < t)
        {
            t = d;
            id = i;
        }
    }
    return {t, id, t < inf};
}
Vec radiance(const Ray& r, int depth, std::mt19937& gen)
{
    std::uniform_real_distribution<double> dis(0.0F, 1.0F);
    //double t = 1e20; // distance to intersection
    //int id = 0;      // id of intersected object
    const Intersection intersection = intersect(r);
    if (!intersection.hit)
    {
        return Vec(); // if miss, return black
    }
    const Sphere& obj = spheres[static_cast<uint32_t>(intersection.id)]; // the hit object
    const Vec x = r.o + r.d * intersection.t;
    const Vec n = (x - obj.p).norm();
    const Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    Vec f = obj.c;
    const double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
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
        const double r1 = 2 * M_PI * dis(gen);
        const double r2 = dis(gen);
        const double r2s = sqrt(r2);
        const Vec w = nl;
        const Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
        const Vec v = w % u;
        const Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, gen));
    }
    if (obj.refl == SPEC) // Ideal SPECULAR reflection
    {
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, gen));
    }
    const Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    const bool into = n.dot(nl) > 0;                // Ray from outside going in?
    constexpr double nc = 1;
    constexpr double nt = 1.5;
    const double nnt = into ? nc / nt : nt / nc;
    const double ddn = r.d.dot(nl);
    const double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    if (cos2t < 0) // Total internal reflection
    {
        return obj.e + f.mult(radiance(reflRay, depth, gen));
    }
    const Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    constexpr double a = nt - nc;
    constexpr double b = nt + nc;
    constexpr double R0 = a * a / (b * b);
    const double c = 1 - (into ? -ddn : tdir.dot(n));
    const double Re = R0 + (1 - R0) * c * c * c * c * c;
    const double Tr = 1 - Re;
    const double P = .25 + .5 * Re;
    const double RP = Re / P;
    const double TP = Tr / (1 - P);
    return obj.e + f.mult(depth > 2 ? (dis(gen) < P ? // Russian roulette
                                           radiance(reflRay, depth, gen) * RP
                                                    : radiance(Ray(x, tdir), depth, gen) * TP)
                                    : radiance(reflRay, depth, gen) * Re + radiance(Ray(x, tdir), depth, gen) * Tr);
}

constexpr Vec AccumulateSampleRadiance(const int samps, const Vec &cx, int sx,
		const uint32_t x, const uint32_t w, const Vec &cy, uint32_t sy,
		const uint32_t y, const uint32_t h, const Ray &cam,
		std::uniform_real_distribution<> &dis, std::mt19937 &gen) {
	Vec r;
	for (int s = 0; s < samps; s++) {
		const double r1 = 2 * dis(gen);
		const double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
		const double r2 = 2 * dis(gen);
		const double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
		const Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5)
				+ cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
		r = r + radiance(Ray(cam.o + d * 140, Vec(d).norm()), 0, gen) * (1. / samps);
	} // Camera rays are pushed ^^^^^ forward to start in interior
	return r;
}

inline constexpr void ProcessPixel(const uint32_t h, const uint32_t y, const uint32_t w,
                         const uint32_t x, const int samps, std::uniform_real_distribution<>& dis,
                         const Vec& cx, const Vec& cy, const Ray& cam, std::mt19937& gen,
                         std::vector<Vec>& c)
{
    for (uint32_t sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
    {
        for (int sx = 0; sx < 2; sx++)
        {
            // 2x2 subpixel cols
            const Vec r = AccumulateSampleRadiance(samps, cx, sx, x, w, cy, sy, y, h, cam, dis, gen);
            c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
        }
    }
}

int main(int argc, char* argv[])
{
    constexpr uint32_t w = 1024;
    constexpr uint32_t h = 768;
    const int samps = argc == 2 ? std::stoi(std::span(argv, static_cast<uint32_t>(argc))[1]) / 4 : 1; // # samples
    constexpr Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());                                  // cam pos, dir
    constexpr Vec cx = Vec(w * .5135 / h);
    constexpr Vec cy = (cx % cam.d).norm() * .5135;
    // Vec r;
    std::vector<Vec> c(static_cast<size_t>(w * h)); // = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1)       // OpenMP
    for (uint32_t y = 0; y < h; y++)
    { // Loop over image rows
        std::stringstream progress;
        progress << "\rRendering (" << samps * 4 << " spp) " << 100. * y / (h - 1);
        std::cerr << progress.str();
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
    for (const auto& printVec : c)
    {
        f << toInt(printVec.x) << " " << toInt(printVec.y) << " " << toInt(printVec.z) << " ";
    }
}
