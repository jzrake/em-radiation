#include <cmath>
#include <utility>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
static const double c = 1.0;

// 3D vector type and operations
struct Vector3D {
    double x;
    double y;
    double z;

    // Array-like access for compatibility
    double& operator[](size_t i) {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range("Vector3D index out of range");
    }

    const double& operator[](size_t i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range("Vector3D index out of range");
    }

    // Calculate vector norm (magnitude)
    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};
// Vector subtraction
static Vector3D operator-(const Vector3D& a, const Vector3D& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// Vector addition
static Vector3D operator+(const Vector3D& a, const Vector3D& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

// Scalar multiplication
static Vector3D operator*(const Vector3D& v, double s) {
    return {v.x * s, v.y * s, v.z * s};
}

// Scalar multiplication (commutative)
static Vector3D operator*(double s, const Vector3D& v) {
    return v * s;
}

// Scalar division
static Vector3D operator/(const Vector3D& v, double s) {
    return {v.x / s, v.y / s, v.z / s};
}

class ParticleTrajectory {
public:
    /** Angular frequency of oscillation */
    double omega = 0.9;

    /** Amplitude of oscillation */
    double ell = 1.0;

    /**
     * Returns the position of a particle at time t
     *
     * The particle oscillates along the z-axis with amplitude ell
     * and angular frequency omega.
     *
     * @param t Time at which to calculate the position
     * @return 3D position vector of the particle
     */
    Vector3D position(double t) const {
        return {0.0, 0.0, ell * std::cos(omega * t)};
    }

    /**
     * Returns the velocity of a particle at time t
     *
     * Calculated as the time derivative of the position function.
     *
     * @param t Time at which to calculate the velocity
     * @return 3D velocity vector of the particle
     */
    Vector3D velocity(double t) const {
        return {0.0, 0.0, -ell * omega * std::sin(omega * t)};
    }
};

// Condition for finding the retarded time
static double retarded_time_condition(const ParticleTrajectory& trajectory, const Vector3D& r, double t, const Vector3D& rprime, double tprime) {
    return (t - tprime) - (r - rprime).norm() / c;
}

// Function to find the retarded time using numerical root finding
static double find_retarded_time(const ParticleTrajectory& trajectory, const Vector3D& r, double t) {
    // Initial guesses for retarded time
    double t0 = t - r.norm() / c;
    double t1 = t;

    // Calculate function values at initial guesses
    double f0 = retarded_time_condition(trajectory, r, t, trajectory.position(t0), t0);
    double f1 = retarded_time_condition(trajectory, r, t, trajectory.position(t1), t1);

    // If initial guesses don't bracket the root, adjust t0
    while (f0 * f1 > 0) {
        t0 = t - trajectory.ell;
        f0 = retarded_time_condition(trajectory, r, t, trajectory.position(t0), t0);
    }

    // Secant method algorithm
    const double tolerance = 1e-10;
    const int max_iterations = 100;

    for (int i = 0; i < max_iterations; i++) {
        // Secant method formula: t_{n+1} = t_n - f(t_n) * (t_n - t_{n-1}) / (f(t_n) - f(t_{n-1}))
        double t2 = t1 - f1 * (t1 - t0) / (f1 - f0);
        double f2 = retarded_time_condition(trajectory, r, t, trajectory.position(t2), t2);

        // Check for convergence
        if (std::abs(f2) < tolerance || std::abs(t2 - t1) < tolerance) {
            return t2;
        }

        // Update for next iteration
        t0 = t1;
        f0 = f1;
        t1 = t2;
        f1 = f2;
    }

    // Return best estimate if max iterations reached
    return t1;
}

/**
 * Calculate the Liénard-Wiechert potentials for a moving charge
 *
 * The Liénard-Wiechert potentials describe the electromagnetic field of a moving point charge.
 * They account for the retardation effects due to the finite speed of light propagation.
 *
 * The scalar potential Φ and vector potential A are calculated as:
 *   Φ = 1/(4πε₀) * q/[|R|*(1-n·v/c)]
 *   A = (v/c) * Φ
 *
 * where:
 *   - R is the vector from the charge's retarded position to the field point
 *   - n is the unit vector in the direction of R
 *   - v is the charge's velocity at the retarded time
 *   - c is the speed of light
 *
 * @param trajectory The trajectory object describing the particle's motion
 * @param r The position vector where the potentials are evaluated
 * @param t The time at which the potentials are evaluated
 * @return A pair containing the scalar potential (first) and vector potential (second)
 */
static std::pair<double, Vector3D> lienard_wiechert_potentials(const ParticleTrajectory& trajectory, const Vector3D& r, double t) {
    // Find the retarded time
    double t_ret = find_retarded_time(trajectory, r, t);

    // Get particle position at retarded time
    Vector3D r_ret = trajectory.position(t_ret);

    // Calculate R = r - r'(t_ret)
    Vector3D R = r - r_ret;
    double R_mag = R.norm();

    // Calculate velocity at retarded time: v = d(r')/dt at t_ret
    Vector3D v = trajectory.velocity(t_ret);

    // Calculate n = R/|R|
    Vector3D n = R / R.norm();

    // Calculate 1 - n·v/c
    double beta_dot_n = (n.x*v.x + n.y*v.y + n.z*v.z) / c;
    double kappa = 1.0 - beta_dot_n;

    // Scalar potential Φ
    double scalar_potential = 1.0 / (4.0 * M_PI * R_mag * kappa);

    // Vector potential A
    Vector3D vector_potential = v * scalar_potential / c;
    return {scalar_potential, vector_potential};
}

PYBIND11_MODULE(em_radiation, m) {
    m.doc() = "Python extension module for calculating radiation from a point charge";

    py::class_<Vector3D>(m, "Vector3D", py::buffer_protocol())
        .def(py::init<>([](double x, double y, double z) {
            return Vector3D{x, y, z};
        }))
        .def("__getitem__", [](const Vector3D &v, size_t i) {
            if (i >= 3) throw py::index_error();
            return v[i];
        })
        .def("__setitem__", [](Vector3D &v, size_t i, double value) {
            if (i >= 3) throw py::index_error();
            v[i] = value;
        })
        .def("__len__", [](const Vector3D&) { return 3; })
        .def("__repr__", [](const Vector3D &v) {
            return "Vector3D(" + std::to_string(v[0]) + ", " +
                                 std::to_string(v[1]) + ", " +
                                 std::to_string(v[2]) + ")";
        })
        .def_property_readonly("norm", &Vector3D::norm, "Calculate vector norm (magnitude)");

    py::class_<ParticleTrajectory>(m, "ParticleTrajectory")
        .def(py::init([](py::kwargs kwargs) {
            ParticleTrajectory trajectory;
            if (kwargs.contains("omega")) trajectory.omega = kwargs["omega"].cast<double>();
            if (kwargs.contains("ell")) trajectory.ell = kwargs["ell"].cast<double>();
            return trajectory;
        }))
        .def_readwrite("omega", &ParticleTrajectory::omega)
        .def_readwrite("ell", &ParticleTrajectory::ell)
        .def("position", &ParticleTrajectory::position, "Returns the position of a particle at time t",
             py::arg("t"))
        .def("velocity", &ParticleTrajectory::velocity, "Returns the velocity of a particle at time t",
             py::arg("t"));

    m.def("retarded_time_condition", &retarded_time_condition,
          "Condition for finding the retarded time",
          py::arg("trajectory"), py::arg("r"), py::arg("t"), py::arg("rprime"), py::arg("tprime"));

    m.def("find_retarded_time", &find_retarded_time,
          "Find the retarded time using numerical root finding",
          py::arg("trajectory"), py::arg("r"), py::arg("t"));

    m.def("lienard_wiechert_potentials", &lienard_wiechert_potentials,
          "Calculate the Liénard-Wiechert potentials for a moving charge",
          py::arg("trajectory"), py::arg("r"), py::arg("t"));

    auto f = [] {};
}
