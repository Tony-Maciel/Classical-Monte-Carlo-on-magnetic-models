/* COMPILE WITH:
   g++ -fopenmp -O3 -march=native -std=c++17 NNIH_correlated.cpp -o NNIH

   RUN WITH:
   ./NNIH [int]
*/

/**************************************************************************************************************************************
 * This program simulates the classical 3D Ising model on the simple hexagonal lattice (honeycomb planes stacked along z, with PBC)
 * with NN, 2NN and 3NN in-plane couplings plus a perpendicular interplane coupling. The Hamiltonian is:
 *
 * H = -\sum_{i,j} J_{i,j} S_i S_j ,  S_i = +/- 1 ,  J_NN = J1, J_2NN = J2, J_3NN = J3, J_perpendicular = Jp.
 *
 * With J1 = +1, J2 = -2, J3 = -2, Jp = +1 the ground state is the 'zig-zag' phase (see Wildes et al., PRB 101, 024415 (2020),
 * Fig. 2 and Table II). The zig-zag state has three symmetry-related domains. Thus, a three-component order parameter O is used.
 *
 * It has Z3 x Z2 symmetry (Lattice + Ising spin). It's expected that a quenched bond disorder creates an implectic phase out of the
 * zig-zag phase in the clean system. This entails spins forming in a zig-zag fashion (3 possible directions) that percolate the
 * lattice ("implexus") and which later break the remaining Z2 spin symmetry. So the Ising symmetry is broken but the lattice symmetry 
 * remains unbroken at low temperatures.
 *
 * The transition temperature for the clean model is around (k_B T_c / J_1) = 8.5 
 *
 * A cold start (random_start = 0 in the input file) prepares a perfect zig-zag ground state. The hot start is randomly oriented spins.
 *
 * --------------------------------------|
 * BOND DISORDER (NEMATIC RANDOM FIELDS) |
 * --------------------------------------|
 * Disorder is implemented as a quenched random bond disorder that couples to the nematic order parameter. This is done by adding 
 * the following term to the Hamiltonian: 
 *
 *  \mathcal{H}_{RF} = -\sum_{i, j, S, a} \phi_{a}(i, j, S) \eta_{a}(i, j, S) , 
 *
 *  the sum over i,j,S is over all sites on the honeycomb lattice (S represents the 2 Bravais lattice base sites: A or B) and the
 *  sum over "a" is over the 3 components of \phi and \eta. Here \eta is the local nematic order parameter given by: 
 * (S_{A,B}(i,j) is the Ising spin value at site (i,j) in that honeycomb plane)
 *
 *  \eta_{1}(i, j, A) = S_A(i,j) * (S_B(i-1,j) + S_B(i,j))
 *  \eta_{1}(i, j, B) = S_B(i,j) * (S_A(i,j)   + S_A(i+1,j))
 *  \eta_{2}(i, j, A) = S_A(i,j) * (S_B(i,j-1) + S_B(i,j))
 *  \eta_{2}(i, j, B) = S_B(i,j) * (S_A(i,j)   + S_A(i,j+1))
 *  \eta_{3}(i, j, A) = S_A(i,j) * (S_B(i,j-1) + S_B(i-1,j))
 *  \eta_{3}(i, j, B) = S_B(i,j) * (S_A(i+1,j) + S_A(i,j+1))
 *
 * While the \phi are random variables drawn from a distribution that is symmetric under permutation. The distribution chosen is: 
 *
 * P(\phi_1, \phi_2, \phi_3) = 1/3 * [\delta(\phi_1 - W)*\delta(\phi_2 + W/2)*\delta(\phi_3 + W/2)
 *                                    +\delta(\phi_1 + W/2)*\delta(\phi_2 - W)*\delta(\phi_3 + W/2)
 *                                    +\delta(\phi_1 + W/2)*\delta(\phi_2 + W/2)*\delta(\phi_3 - W)]
 *
 * Basically, if you expand out the random field Hamiltonian by using the above definition of the local nematic order parameter and group 
 * up terms in the original Hamiltonian, you'll find that the result of this disorder is just to change the values of the J1 bonds. This 
 * change essentially weakens one of the J1 bonds and strengthens the other 2 J1 bonds, favoring a certain zig zag direction for that site.
 *
 * To correlate these random nematic fields, first define a kernel size (up to NN distances, up to NNN distances etc.). Then, for each site, 
 * generate 2 independent Gaussian random numbers (mean = 0, variance = 1) which represent the real and imaginary parts of a complex number.
 * Then, go to each site, sum these complex Gaussian numbers over the kernel size centered at that site, divide by the square root of the total 
 * number of sites considered (normalization), calculate the argument of the resulting complex number and choose one of the 3 possible permutations
 * for \phi given above depending on the angle this argument makes (e.g., between 0 and 2\pi/3, choose \phi_1=W, \phi_2=-W/2 and \phi_3=-W/2, if 
 * between 2\pi/3 and 4\pi/3, choose ...).  
 *
 *   kernel | sites | layer k             | layers k+-1          | layers k+-2 | R
 *   -------+-------+---------------------+----------------------+-------------+----------
 *      0   |   1   | (uncorrelated)      |                      |             |
 *      1   |   6   | site + 1NN          | site above/below     |             | a
 *      2   |  10   | site + 1NN + 2NN    |                      |             | in-plane
 *      3   |  13   | site + 1NN..3NN     |                      |             | in-plane
 *      4   |  18   | site + 1NN + 2NN    | site + its 1NN       |             | sqrt(3) a
 *      5   |  35   | site + 1NN..3NN     | site + its 1NN + 2NN | site        | 2a
 *
 * Kernels 1, 4 and 5 are the full spheres of radius a, sqrt(3)a and 2a (boundary included). Note that R is
 * the radius of the summed cluster: two sites share part of their clusters, and hence have correlated
 * disorder, out to 2R.
 *
 * The goal of correlating the disorder is to give more room for the implectic order to appear.
 *
 *-----------------| 
 * IMPLECTIC ORDER |
 * ----------------|
 * At low temperatures, it's expected that a phase appears that breaks the Z2 Ising spin symmetry, but keeps the Z3 lattice symmetry unbroken.
 * An interwoven ("implexus") network of the 3 possible zig zag ground states percolate the lattice, and each of these percolating clusters 
 * individually breaks the remaining Z2 Ising spin symmetry.
 *
 * -----------------|
 * TEMPERATURE GRID |
 * -----------------|
 * Every system size uses exactly n_rep temperatures on [Ti, Tf] (both ends included, evenly spaced in beta),
 * so all L are simulated at the same temperature points. This is done so that when calculating the disorder 
 * averages, every temperature point has the same number of disorder samples. The beta spacing needed for good
 * swap acceptance shrinks roughly as L^(-3/2) (1.5 = (dimension of system) / 2 , as suggested by Hukushima and Nemoto.),
 * so n_rep has to be chosen for the LARGEST L.
 *
 * Each thread handles length/Nthreads replicas, so length must be a multiple of Nthreads. When the count on
 * [Ti, Tf] is not, the missing replicas are added ABOVE Tf, continuing the same beta spacing, and the points
 * on [Ti, Tf] are left exactly where they were. Those padding temperatures take part in parallel tempering
 * and appear in <o>_pt.dat and <o>_ehist.dat, but <o>.dat only has the rows on [Ti, Tf], so it has the
 * same temperature column for every L and every thread count. 
 *
 * N.B.: In my experience, it's ideal to: not use too many replicas; have around 40% swap acceptance rate; swap often; 
 * extend a considerable amount to the paramagnetic phase to help with equilibration issues at lower temperatures.
 *
 * ----------|
 * SNAPSHOTS |
 * ----------|
 * With snap_every > 0 in the input file, the spin configuration of whichever replica currently sits at
 * the LOWEST temperature is appended to <o>_snap.dat every snap_every parallel-tempering rounds.
 *************************************************************************************************************************************/

#include <iostream>   // Input and output
#include <vector>     // Dynamic arrays
#include <array>      // static arrays
#include <cmath>      // exp
#include <fstream>    // Reading and writing to files
#include <string>     // Unique file names
#include <iomanip>    // Set precision of output data
#include <algorithm>  // Std::fill (a little faster than for loop)
#include <cstdint>    // Fixed width integer types (in prng)
#include <cstddef>    // std::size_t
#include <limits>     // numeric_limits (energy histogram range)
#include <memory>     // std::unique_ptr (the System object is far too large for the stack)
#include <chrono>     // Getting current time
#include <ctime>      // Making time from chrono more human readable
#include <sstream>    // Treats strings as files
#include <omp.h>      // Parallelize
#include <filesystem> // To check if file exists

using RealType = double;
using IntType = long long;

uint64_t master_state[4];

/* xoshiro256+ implementation by David Blackman and Sebastiano Vigna (vigna@acm.org)
 * Public domain.
 *
 * The period of this prng is 2^256 - 1 (roughly 10^77). jump() is being used for each thread
 * at a unique temperature, then long_jump() is used at a certain realization of disorder
 * to guarantee new random numbers.
 */
struct Xoshiro256Plus {
    uint64_t state[4];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    // splitmix64 for seeding as recommended by the authors of xoshiro256+
    void seed(uint64_t seed_val) {
        uint64_t z = seed_val;
        for (int i = 0; i < 4; i++) {
            z += 0x9e3779b97f4a7c15;
            uint64_t s = z;
            s = (s ^ (s >> 30)) * 0xbf58476d1ce4e5b9;
            s = (s ^ (s >> 27)) * 0x94d049bb133111eb;
            state[i] = s ^ (s >> 31);
        }
    }

    uint64_t next() {
        const uint64_t result = state[0] + state[3];
        const uint64_t t = state[1] << 17;
        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];
        state[2] ^= t;
        state[3] = rotl(state[3], 45);
        return result;
    }

    double next_double() {
        return (double)(next() >> 11) * 0x1.0p-53;
    }

    /* Two independent standard normals per call (Box-Muller). Used only when building spatially
     * correlated disorder.
     */
    void gaussian_pair(double& z1, double& z2) {
        double u1 = next_double();
        if (u1 < 1e-300) // safety check because next_double() can, in principle, return 0
            u1 = 1e-300;
        const double u2 = next_double();
        const double r  = std::sqrt(-2.0 * std::log(u1));
        const double th = 6.283185307179586476925 * u2;
        z1 = r * std::cos(th);
        z2 = r * std::sin(th);
    }

    // Advances 2^128 steps (roughly 10^38). Used for threads within a temperature range.
    void jump() {
        static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < 4; i++) {
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & 1ULL << b) {
                    s0 ^= state[0];
                    s1 ^= state[1];
                    s2 ^= state[2];
                    s3 ^= state[3];
                }
                next();
            }
        }
        state[0] = s0;
        state[1] = s1;
        state[2] = s2;
        state[3] = s3;
    }

    // Advances 2^192 steps (roughly 10^57). Used to move to the next realization of disorder safely.
    void long_jump() {
        static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe69d };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < 4; i++) {
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & 1ULL << b) {
                    s0 ^= state[0];
                    s1 ^= state[1];
                    s2 ^= state[2];
                    s3 ^= state[3];
                }
                next();
            }
        }
        state[0] = s0;
        state[1] = s1;
        state[2] = s2;
        state[3] = s3;
    }
};

inline int fast_mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

/* Honeycomb-in-plane geometry, expressed as offsets in the triangular Bravais cell indices (i,j).
 * Cell (i,j) holds two basis sites: s = 0 (sublattice A) at i*a1 + j*a2, and s = 1 (sublattice B)
 * at the same point + delta, with a1 = (3/2, sqrt(3)/2), a2 = (3/2, -sqrt(3)/2), delta = (1,0).
 *
 *   1NN : 3 sites at distance a,        on the OPPOSITE sublattice
 *   2NN : 6 sites at distance sqrt(3)a, on the SAME sublattice
 *   3NN : 3 sites at distance 2a,       on the OPPOSITE sublattice
 *
 * The tables below are written for sublattice A. For sublattice B every offset is negated
 * (o2 is inversion symmetric, so it is identical for both).
 *
 *   slot 0 : A(i,j) - B(i-1,j)   frustrated in psi_1  [ psi_1 = sum (-1)^i     (S_A + S_B) ]
 *   slot 1 : A(i,j) - B(i,j-1)   frustrated in psi_2  [ psi_2 = sum (-1)^j     (S_A + S_B) ]
 *   slot 2 : A(i,j) - B(i,j)     frustrated in psi_3  [ psi_3 = sum (-1)^(i+j) (S_A - S_B) ]
 */
inline constexpr int o1[3][2] = { {-1, 0}, { 0,-1}, { 0, 0} };
inline constexpr int o2[6][2] = { { 1, 0}, {-1, 0}, { 0, 1}, { 0,-1}, { 1,-1}, {-1, 1} };
inline constexpr int o3[3][2] = { {-1,-1}, {-1, 1}, { 1,-1} };

// represents an instance of the system being simulated
template <int L>
struct alignas(64) System {

    // L cells along a1, L cells along a2, L layers along z; each cell carries TWO sites (honeycomb basis).
    static constexpr int Ncell = L * L * L;
    static constexpr int Total = 2 * Ncell;

    /* Coordination number Z, laid out contiguously in Neighbors:
     *   [ 0,  3)  1NN (3 sites)           -> J1
     *   [ 3,  9)  2NN (6 sites)           -> J2
     *   [ 9, 12)  3NN (3 sites)           -> J3
     *   [12, 14)  perpendicular (2 sites) -> Jp
     */
    static constexpr int Z = 14;

    // size_t instead of int to avoid overflow in some template cases
    static constexpr std::size_t NbSize = std::size_t(Total) * std::size_t(Z);

    static constexpr int NExpA = 216;
    static constexpr int NExpB = 84;

    inline static std::array<int, NbSize> Neighbors{};  // address of every neighbor of every site (shared)
    std::array<int8_t, Total> Status;                   // Ising spin (+1 or -1) on each site (private to each replica)
    inline static std::array<RealType, Z> Jbond{};      // coupling attached to each neighbor slot

    /* Jnn[site*3 + m] is the class of the bond in 1NN slot m at that site, whose
     * coupling is Jval[Jnn[site*3 + m]].
     *
     * Both endpoints of a bond store the same number. In the clean model every Jval entry equals J1;
     * with nematic random fields the entries become bond dependent.
     *
     * N.B.: Neighbors and Jnn are IDENTICAL for every replica (one lattice and one disorder
     * realization per Run()). They are therefore SHARED by all System<L> instances (C++17 inline static),
     * written once before the parallel region and read-only from then on to avoid race conditions.
     *
     * Because of this, init_sys(), set_neighbors() and set_nematic_disorder() must be called on ONE
     * System<L> only, before the replicas are created. Calling them on a replica would silently
     * rewrite the lattice underneath every other replica.
     */
    inline static std::array<uint8_t, std::size_t(Total) * 3> Jnn{};
    inline static std::array<RealType, 3> Jval{};

    inline static std::vector<RealType> ExpA{};
    inline static std::vector<RealType> ExpB{};

    Xoshiro256Plus rng;

    // auxiliary functions for readability...
    inline int spin(int i) const {
        return int(Status[i]);
    }

    inline void flip(int i) {
        Status[i] = int8_t(-Status[i]);
    }

    inline void put(int i, bool down) {
        Status[i] = down ? int8_t(-1) : int8_t(1);
    }

    void init_sys(bool random_start, RealType J1, RealType J2, RealType J3, RealType Jp,
                  int zigzag_domain = 0) {
        for (int m = 0;  m < 3;  ++m) // NN 
            Jbond[m] = J1;
        for (int m = 3;  m < 9;  ++m) // second NN 
            Jbond[m] = J2;
        for (int m = 9;  m < 12; ++m) // third NN 
            Jbond[m] = J3;
        for (int m = 12; m < 14; ++m) // up and down neighbors (among planes)
            Jbond[m] = Jp;

        // clean value everywhere
        // set_nematic_disorder() overwrites this when W > 0
        std::fill(Jnn.begin(), Jnn.end(), uint8_t(0));
        std::fill(Jval.begin(), Jval.end(), J1);

        if (!random_start) {
            set_zigzag(zigzag_domain); // cold start: a perfect zig-zag ground state
        } else {
            for (int i = 0; i < Total; ++i)
                put(i, rng.next_double() < 0.5); // hot start
        }
    }

    /* Fills the lattice with one of the three zig-zag ground states.
     * The domain argument selects which of the three symmetry related states is used, and matches
     * the order parameter components defined in Get():
     *
     *   domain 0 : S_A(i,j,k) = S_B(i,j,k) =  (-1)^i        -> psi_1 = 1, the 1NN slot-0 bonds are frustrated
     *   domain 1 : S_A(i,j,k) = S_B(i,j,k) =  (-1)^j        -> psi_2 = 1, the 1NN slot-1 bonds are frustrated
     *   domain 2 : S_A(i,j,k) = -S_B(i,j,k) = (-1)^(i+j)    -> psi_3 = 1, the 1NN slot-2 bonds are frustrated
     *
     * The state is uniform along z, which is what Jp > 0 wants.
     * (The (-1)^i and (-1)^j patterns only close on themselves if L is even.)
     */
    void set_zigzag(int domain) {
        static_assert(L % 2 == 0, "The zig-zag cold start requires an even L to match the periodic boundaries.");

        if (domain < 0 || domain > 2)
            domain = 0;

        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    int parity;
                    if (domain == 0)
                        parity = i % 2;
                    else if (domain == 1)
                        parity = j % 2;
                    else
                        parity = (i + j) % 2;

                    const int8_t S_A = parity ? int8_t(-1) : int8_t(1);
                    const int8_t S_B = (domain == 2) ? int8_t(-S_A) : S_A;

                    Status[idx(i, j, k, 0)] = S_A;
                    Status[idx(i, j, k, 1)] = S_B;
                }
            }
        }
    }

    /* Maps a lattice site to a positive integer.
     * (i, j) are the triangular Bravais cell indices, k is the layer, s is the sublattice (A=0, B=1).
     * Periodic boundary conditions are applied here.
     */
    static inline int idx(int i, int j, int k, int s) {
        i = fast_mod(i, L);
        j = fast_mod(j, L);
        k = fast_mod(k, L);
        return 2 * ((k * L + i) * L + j) + s;
    }

    /* 
     * Builds the neighbor table for the simple hexagonal (stacked honeycomb) lattice.
     */
    void set_neighbors() {
        for (int k = 0; k < L; ++k) { // for each z plane 
            for (int i = 0; i < L; ++i) { // for each row in the honeycomb lattice 
                for (int j = 0; j < L; ++j) { // for each column in the honeycomb lattice 
                    for (int s = 0; s < 2; ++s) { // for each basis point in the triangular Bravais lattice
                        const int site = idx(i, j, k, s);
                        const int g    = s ? -1 : 1; // sublattice B has s=1, so g=1, while A has s=0 and g=-1
                        int* nb        = &Neighbors[std::size_t(site) * Z];

                        // NN: 3 sites, opposite sublattice (1-s)
                        for (int m = 0; m < 3; ++m)
                            nb[m] = idx(i + g*o1[m][0], j + g*o1[m][1], k, 1 - s);

                        // second NN: 6 sites, same sublattice
                        for (int m = 0; m < 6; ++m)
                            nb[3 + m] = idx(i + o2[m][0], j + o2[m][1], k, s);

                        // third NN: 3 sites, opposite sublattice (1-s)
                        for (int m = 0; m < 3; ++m)
                            nb[9 + m] = idx(i + g*o3[m][0], j + g*o3[m][1], k, 1 - s);

                        // Perpendicular: 2 sites, same sublattice
                        nb[12] = idx(i, j, k + 1, s);
                        nb[13] = idx(i, j, k - 1, s);
                    }
                }
            }
        }
    }

    /* Used only to build the 3D correlation kernels (4 and 5) of the nematic disorder.
     *
     * Appends to 'out' the sites of layer k+dk that lie in the in-plane shells 0..nmax around the
     * column of site (i, j, k, s):
     *   shell 0 : the site itself (or, for dk != 0, the site directly above/below it)   1 site
     *   shell 1 : 1NN, distance a,        opposite sublattice                             3 sites
     *   shell 2 : 2NN, distance sqrt(3)a, same sublattice                                 6 sites
     *   shell 3 : 3NN, distance 2a,       opposite sublattice                             3 sites
     * (The offsets and the sign g are exactly those of set_neighbors().)
     */
    static void append_shells(int i, int j, int k, int s, int dk, int nmax, std::vector<int>& out) {
        const int g  = s ? -1 : 1;
        const int kk = k + dk;

        if (nmax >= 0)
            out.push_back(idx(i, j, kk, s));
        if (nmax >= 1)
            for (int m = 0; m < 3; ++m)
                out.push_back(idx(i + g*o1[m][0], j + g*o1[m][1], kk, 1 - s));
        if (nmax >= 2)
            for (int m = 0; m < 6; ++m)
                out.push_back(idx(i + o2[m][0], j + o2[m][1], kk, s));
        if (nmax >= 3)
            for (int m = 0; m < 3; ++m)
                out.push_back(idx(i + g*o3[m][0], j + g*o3[m][1], kk, 1 - s));
    }

    /* Outermost in-plane shell kept in layer k+dz, indexed by |dz| = 0, 1, 2 (-1 = layer not used).
     * With the interlayer spacing equal to a, a site in layer k+dz at in-plane distance d from the column
     * is at distance sqrt(d^2 + dz^2) a, so these are the spheres of radius sqrt(3)a and 2a:
     *
     *   kernel 4 : |dz|=0 -> 2NN (sqrt3 a), |dz|=1 -> 1NN (sqrt2 a)                   1+3+6   + 2*(1+3)            = 18
     *   kernel 5 : |dz|=0 -> 3NN (2a),      |dz|=1 -> 2NN (2a),  |dz|=2 -> site (2a)  1+3+6+3 + 2*(1+3+6) + 2*1    = 35
     */
    static constexpr int kernel4_shells[3] = { 2, 1, -1 };
    static constexpr int kernel5_shells[3] = { 3, 2,  0 };

    // Every site of the 3D kernel centred on (i, j, k, s): its own layer first, then k+1, k-1, k+2, k-2
    static void kernel_cluster(int i, int j, int k, int s, const int* shells, std::vector<int>& out) {
        out.clear();
        for (int adz = 0; adz <= 2; ++adz) {
            if (shells[adz] < 0)
                continue;
            append_shells(i, j, k, s, +adz, shells[adz], out);
            if (adz > 0)
                append_shells(i, j, k, s, -adz, shells[adz], out);
        }
    }

    /* Builds one quenched realization of the nematic random field.
     * must be called AFTER set_neighbors(), and BEFORE the replicas are copied out of this object.
     */
    void set_nematic_disorder(RealType J1, RealType W, int kernel) {
        if (W == 0.0) return; // clean model: Jval already holds J1 everywhere

        // favored zig-zag direction at every site; phi is then W on that slot and -W/2 on the other two
        const std::size_t ntot = std::size_t(Total);
        std::vector<int8_t> dir(ntot);

        if (kernel <= 0) {
            /* uncorrelated disorder: one independent draw per site
             */
            for (std::size_t p = 0; p < ntot; ++p)
                dir[p] = int8_t(rng.next_double() * 3.0); // choose one permutation
        } else {
            /* correlated disorder
             * Neighbors matrix is laid out contiguously, so "site + everything out to neighbor n" is
             * simply slots [0, nslots) plus the site itself.
             */
            std::vector<RealType> gR(ntot), gI(ntot);

            // Generate one complex Gaussian distributed number per site
            for (std::size_t p = 0; p < ntot; ++p) {
                RealType z1, z2;
                rng.gaussian_pair(z1, z2);
                gR[p] = z1;
                gI[p] = z2;
            }

            std::vector<RealType> hR(ntot), hI(ntot);

            if (kernel <= 3) {
                const int      nslots     = (kernel == 1) ? 3 : (kernel == 2) ? 9 : 12;
                const bool     use_perp   = (kernel == 1);
                const int      nterms     = nslots + 1 + (use_perp ? 2 : 0); // use_perp is used because I thought that the correlation was only in plane...
                const RealType inv_sqrt_n = 1.0 / std::sqrt(RealType(nterms));

                for (int site = 0; site < Total; ++site) {
                    // get complex Gaussian number at that site...
                    const int* nb = &Neighbors[std::size_t(site) * Z];
                    RealType sR = gR[std::size_t(site)];
                    RealType sI = gI[std::size_t(site)];

                    // ...and sum it with neighbors
                    for (int m = 0; m < nslots; ++m) {
                        sR += gR[std::size_t(nb[m])];
                        sI += gI[std::size_t(nb[m])];
                    }

                    if (use_perp) { // get perpendicular neighbors too...
                        sR += gR[std::size_t(nb[12])] + gR[std::size_t(nb[13])];
                        sI += gI[std::size_t(nb[12])] + gI[std::size_t(nb[13])];
                    }

                    // normalize
                    hR[std::size_t(site)] = inv_sqrt_n * sR;
                    hI[std::size_t(site)] = inv_sqrt_n * sI;
                }
            } else {
                /* Kernels 4 (18 sites, R = sqrt(3)a) and 5 (35 sites, R = 2a). 
                 * note that the clusters need sites (1NN and 2NN of the sites above/below,
                 * sites two layers away) that are not in Neighbors, so they are built on the fly by kernel_cluster().
                 */
                const int* shells = (kernel == 4) ? kernel4_shells : kernel5_shells;
                std::vector<int> members;
                members.reserve(64);

                // A cluster must not wrap onto itself through the periodic boundaries (kernel 5 needs L >= 5)
                bool wraps = false;
                for (int s = 0; s < 2; ++s) {
                    kernel_cluster(0, 0, 0, s, shells, members);
                    std::sort(members.begin(), members.end());
                    if (std::adjacent_find(members.begin(), members.end()) != members.end())
                        wraps = true;
                }
                if (wraps)
                    std::cerr << "Warning: L=" << L << " is too small for kernel " << kernel
                              << ": the cluster wraps around the periodic boundaries, so some sites are summed twice."
                              << std::endl;

                const RealType inv_sqrt_n = 1.0 / std::sqrt(RealType(members.size()));

                for (int k = 0; k < L; ++k) {
                    for (int i = 0; i < L; ++i) {
                        for (int j = 0; j < L; ++j) {
                            for (int s = 0; s < 2; ++s) {
                                kernel_cluster(i, j, k, s, shells, members);

                                RealType sR = 0.0, sI = 0.0;
                                for (const int q : members) {
                                    sR += gR[std::size_t(q)];
                                    sI += gI[std::size_t(q)];
                                }

                                const std::size_t site = std::size_t(idx(i, j, k, s));
                                hR[site] = inv_sqrt_n * sR;
                                hI[site] = inv_sqrt_n * sI;
                            }
                        }
                    }
                }
            }

            /*
             * Now get the argument of the resulting complex number (sum of complex gaussians) and map 
             * it to one of the three 120 degree arcs on the unit circle on the complex plane. Each sector 
             * corresponds to one chosen permutation.
             * */
            constexpr RealType two_pi = 6.283185307179586476925; // a slice of pi... 

            for (int site = 0; site < Total; ++site) {
                RealType th = std::atan2(hI[std::size_t(site)], hR[std::size_t(site)]);
                if (th < 0.0) // wrap around to unit circle
                    th += two_pi;

                int d = int(th * (3.0 / two_pi)); // get corresponding 120 degree sector
                if (d > 2) // in case of numerical error of floating point numbers
                    d = 2;          
                dir[std::size_t(site)] = int8_t(d);
            }
        }

        // phi_m at a site: W on the favored slot, -W/2 on the other two
        Jval[0] = J1 + W;        // J1 -> J1 + W + W - W/2 - W/2 
        Jval[1] = J1 - 0.5 * W;  // J1 -> J1 + W - W/2 - W/2 - W/2
        Jval[2] = J1 - 2.0 * W;  // J1 -> J1 - W/2 - W/2 - W/2 - W/2

        // Enumerates every NN bond exactly once.
        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    const int siteA = idx(i, j, k, 0);

                    for (int m = 0; m < 3; ++m) {
                        const int siteB = Neighbors[std::size_t(siteA) * Z + m]; // partner, same slot m

                        const uint8_t J_bond = uint8_t(int(dir[std::size_t(siteA)] == m)
                                                     + int(dir[std::size_t(siteB)] == m));

                        Jnn[std::size_t(siteA) * 3 + m] = J_bond;
                        Jnn[std::size_t(siteB) * 3 + m] = J_bond;
                    }
                }
            }
        }
    }

    // Local exchange field h_i = sum_j J_ij S_j, summed over all Z neighbors of the site
    inline RealType local_field(int site) const {
        const int* nb = &Neighbors[std::size_t(site) * Z];
        RealType h = 0.0;

        // 1NN: bond resolved coupling (carries the nematic random field)
        for (int nn = 0; nn < 3; ++nn)
            h += Jval[Jnn[std::size_t(site) * 3 + nn]] * RealType(spin(nb[nn]));

        // 2NN, 3NN and perpendicular
        for (int nn = 3; nn < Z; ++nn)
            h += Jbond[nn] * RealType(spin(nb[nn]));

        return h;
    }

    // beta[k] is the inverse temperature of slot k (Measurements::beta_of)
    static void build_exp_tables(const std::vector<RealType>& beta) {
        const int length = int(beta.size());
        ExpA.assign(std::size_t(length) * NExpA, 0.0);
        ExpB.assign(std::size_t(length) * NExpB, 0.0);

        for (int k = 0; k < length; ++k) {
            const RealType mb2 = -2.0 * beta[std::size_t(k)];

            for (int c2 = 0; c2 < 6; ++c2) {
                for (int c1 = 0; c1 < 6; ++c1) {
                    for (int c0 = 0; c0 < 6; ++c0) {
                        const RealType v = (c0 < 3 ? Jval[c0] : -Jval[c0 - 3])
                                         + (c1 < 3 ? Jval[c1] : -Jval[c1 - 3])
                                         + (c2 < 3 ? Jval[c2] : -Jval[c2 - 3]);

                        ExpA[std::size_t(k) * NExpA + c0 + 6 * c1 + 36 * c2] =
                            std::exp(std::min(mb2 * v, 700.0));
                    }
                }
            }

            for (int ip = 0; ip < 3; ++ip) {
                for (int i3 = 0; i3 < 4; ++i3) {
                    for (int i2 = 0; i2 < 7; ++i2) {
                        const RealType v = Jbond[3]  * RealType(2 * i2 - 6)
                                         + Jbond[9]  * RealType(2 * i3 - 3)
                                         + Jbond[12] * RealType(2 * ip - 2);

                        ExpB[std::size_t(k) * NExpB + i2 + 7 * (i3 + 4 * ip)] =
                            std::exp(std::min(mb2 * v, 700.0));
                    }
                }
            }
        }
    }

    // One MCS with the Metropolis algorithm
    void MCS_Metropolis(int k_temp) {
        const RealType* eA = &ExpA[std::size_t(k_temp) * NExpA];
        const RealType* eB = &ExpB[std::size_t(k_temp) * NExpB];

        int8_t*        __restrict__ st    = Status.data();
        const int*     __restrict__ nbase = Neighbors.data();
        const uint8_t* __restrict__ jbase = Jnn.data();

        for (int step = 0; step < Total; ++step) {
            int site = int(rng.next_double() * Total);

            const int*     nb = nbase + std::size_t(site) * Z;
            const uint8_t* jc = jbase + std::size_t(site) * 3;
            const int      Si = int(st[site]);

            const int c0 = int(jc[0]) + (Si * int(st[nb[0]]) > 0 ? 0 : 3);
            const int c1 = int(jc[1]) + (Si * int(st[nb[1]]) > 0 ? 0 : 3);
            const int c2 = int(jc[2]) + (Si * int(st[nb[2]]) > 0 ? 0 : 3);

            int r2 = 0;
            for (int nn = 3; nn < 9; ++nn)
                r2 += int(st[nb[nn]]);

            int r3 = 0;
            for (int nn = 9; nn < 12; ++nn)
                r3 += int(st[nb[nn]]);

            const int rp = int(st[nb[12]]) + int(st[nb[13]]);

            // Flipping S_i -> -S_i changes the energy by dE = 2 * S_i * h_i
            const RealType w = eA[c0 + 6 * c1 + 36 * c2]
                             * eB[((Si * r2 + 6) >> 1) + 7 * (((Si * r3 + 3) >> 1) + 4 * ((Si * rp + 2) >> 1))];

            if (w > 1.0 || rng.next_double() < w) {
                st[site] = int8_t(-st[site]);
            }
        }
    }
}
struct Measurements {
    int mcs, pt_mcs, burnin, lag, length, Nthreads;
    int n_rep;               // temperatures on [Ti, Tf], both ends included: the same for every L
    int n_pad;               // extra temperatures ABOVE Tf, only there to make length a multiple of Nthreads
    RealType Bi, Bf, dB, dB_pad, current_beta;
    RealType E1, m2_cs, m4_cs;
    RealType O2_cs, O4_cs;
    RealType psi_norm_cs;    // |psi| = sqrt(psi_1^2 + psi_2^2 + psi_3^2), configuration by configuration
    RealType psi_abs_cs[3];  // |psi_1|, |psi_2|, |psi_3|
    RealType psi_prod_cs;    // |psi_1 * psi_2 * psi_3|
    RealType eta_cs[3];      // the three components of the local nematic order parameter
    RealType Nl2_cs, Nl4_cs; // local (bond-energy) nematic order parameter
    RealType Ng2_cs, Ng4_cs; // global nematic order parameter, built from psi_m^2
    RealType Nl_abs_cs;      // |Nl|, configuration by configuration
    RealType Ng_abs_cs;      // |Ng|, configuration by configuration

    /* Temperature grid
     */
    RealType beta_of(int k) const {
        return (k >= n_pad) ? Bf + dB * RealType(k - n_pad)
                            : Bf - dB_pad * RealType(n_pad - k);
    }

    void init_meas(int mcs_in, int pt_mcs_in, int burnin_in, int lag_in,
                   RealType Ti_in, RealType Tf_in, int n_rep_in) {
        Nthreads = omp_get_max_threads();
        mcs      = mcs_in;
        pt_mcs   = pt_mcs_in;
        burnin   = burnin_in;
        lag      = lag_in;

        n_rep    = n_rep_in;  // main() guarantees n_rep >= 2

        n_pad  = (Nthreads - n_rep % Nthreads) % Nthreads;
        length = n_rep + n_pad;

        Bi = 1.0 / Ti_in;                  // Ti_in is the LOW temperature end -> LARGEST beta
        Bf = 1.0 / Tf_in;                  // Tf_in is the HIGH  temperature end -> SMALLEST  beta
        dB = (Bi - Bf) / static_cast<RealType>(n_rep - 1);
        // temperature slot k has beta = Bf + dB*k, so slot 0 is the hottest and slot length-1 the coldest
        // (with padding above Tf the interval starts at slot n_pad instead: slot k has beta_of(k), see above)

        dB_pad = dB;
        if (n_pad > 0 && dB_pad > Bf / RealType(n_pad + 1))
            dB_pad = Bf / RealType(n_pad + 1);
    }
};

/* Gathers relevant observables.
 *
 * Order parameter components:
 *   psi_1 = (1/N) sum (-1)^i     (S_A + S_B)
 *   psi_2 = (1/N) sum (-1)^j     (S_A + S_B)
 *   psi_3 = (1/N) sum (-1)^(i+j) (S_A - S_B)
 *
 * O^2 = psi_1^2 + psi_2^2 + psi_3^2 equals 1 in any perfect zig-zag ground state.
 *
 * Nematic order parameters: n = n_1 + n_2 * e^(2i pi/3) + n_3 * e^(4i pi/3) .
 * Explicitly, Re(n) = n_1 - n_2/2 - n_3/2 and Im(n) = (sqrt3/2)(n_2 - n_3). 
 *
 * local nematic order parameter Nl: n_m = (1/N) sum_site eta_m(site).
 *
 * global nematic order parameter Ng: n_m = psi_m^2. 
 */
template <int L>
void Get(const System<L>& sys, Measurements& meas) {
    RealType m = 0.0, E1 = 0.0;
    RealType psi1 = 0.0, psi2 = 0.0, psi3 = 0.0;
    RealType eta1 = 0.0, eta2 = 0.0, eta3 = 0.0;

    // first get magnetization, energy and the local nematic operator
    for (int site = 0; site < sys.Total; ++site) {
        RealType S_i = RealType(sys.spin(site));
        m           += S_i;

        // Exchange energy. Every bond is visited twice across the site loop, hence the factor 1/2.
        E1 -= 0.5 * S_i * sys.local_field(site);

        // get the three 1 NN
        const int* nb = &sys.Neighbors[std::size_t(site) * System<L>::Z];
        RealType S_0  = RealType(sys.spin(nb[0]));
        RealType S_1  = RealType(sys.spin(nb[1]));
        RealType S_2  = RealType(sys.spin(nb[2]));

        eta1 += S_i * (S_1 + S_2);
        eta2 += S_i * (S_0 + S_2);
        eta3 += S_i * (S_0 + S_1);
    }

    // now get order parameters
    for (int k = 0; k < L; ++k) { // for every layer... 
        for (int i = 0; i < L; ++i) { // for every row... 
            for (int j = 0; j < L; ++j) { // for every column...
                RealType S_A = RealType(sys.spin(System<L>::idx(i, j, k, 0)));
                RealType S_B = RealType(sys.spin(System<L>::idx(i, j, k, 1)));

                RealType ei  = (i % 2) ? -1.0 : 1.0;
                RealType ej  = (j % 2) ? -1.0 : 1.0;

                psi1 += ei * (S_A + S_B);
                psi2 += ej * (S_A + S_B);
                psi3 += ei * ej * (S_A - S_B);
            }
        }
    }

    RealType one_over_Total = 1.0 / RealType(sys.Total);
    m    *= one_over_Total;
    psi1 *= one_over_Total;
    psi2 *= one_over_Total;
    psi3 *= one_over_Total;
    eta1 *= one_over_Total;
    eta2 *= one_over_Total;
    eta3 *= one_over_Total;

    meas.E1    = E1 * one_over_Total;
    meas.m2_cs = m * m;
    meas.m4_cs = meas.m2_cs * meas.m2_cs;

    meas.psi_abs_cs[0] = std::fabs(psi1);
    meas.psi_abs_cs[1] = std::fabs(psi2);
    meas.psi_abs_cs[2] = std::fabs(psi3);

    meas.psi_prod_cs = std::fabs(psi1 * psi2 * psi3);

    RealType p1 = psi1 * psi1, p2 = psi2 * psi2, p3 = psi3 * psi3;

    meas.O2_cs = p1 + p2 + p3;
    meas.O4_cs = meas.O2_cs * meas.O2_cs;

    meas.psi_norm_cs = std::sqrt(meas.O2_cs);

    constexpr RealType half_sqrt3 = 0.86602540378443864676;

    /* Local (bond energy) nematic. The components are kept with the factor 1/2 folded in, so that
     * n_m = (1/2N) sum_site eta_m(site) lies in [-1, 1] and Nl is exactly
     *     Re(Nl) = n_1 - n_2/2 - n_3/2 ,   Im(Nl) = (sqrt3/2)(n_2 - n_3)
     */
    meas.eta_cs[0] = 0.5 * eta1;
    meas.eta_cs[1] = 0.5 * eta2;
    meas.eta_cs[2] = 0.5 * eta3;

    RealType Nl_re = meas.eta_cs[0] - 0.5 * (meas.eta_cs[1] + meas.eta_cs[2]);
    RealType Nl_im = half_sqrt3 * (meas.eta_cs[1] - meas.eta_cs[2]);
    meas.Nl2_cs    = Nl_re * Nl_re + Nl_im * Nl_im;
    meas.Nl4_cs    = meas.Nl2_cs * meas.Nl2_cs;

    RealType Ng_re = p1 - 0.5 * p2 - 0.5 * p3;
    RealType Ng_im = half_sqrt3 * (p2 - p3);
    meas.Ng2_cs    = Ng_re * Ng_re + Ng_im * Ng_im;
    meas.Ng4_cs    = meas.Ng2_cs * meas.Ng2_cs;

    // For <|Nl|> and <|Ng|> 
    meas.Nl_abs_cs = std::sqrt(meas.Nl2_cs);
    meas.Ng_abs_cs = std::sqrt(meas.Ng2_cs);
}

/*
 * Runs a simulation for a single fixed size (L).
 * */
template <int L>
void Run(int mcs_in, int pt_mcs_in, int burnin_in, int lag_in, RealType Ti_in, RealType Tf_in,
         RealType J1, RealType J2, RealType J3, RealType Jp, RealType W, int kernel,
         int n_rep, int n_ebins, int snap_every, bool random_start, std::string output_file_name) {

    // Heap-allocated: System<L> is far too large to sit on the stack once L grows
    // (the stack is about 8 Mb. Check with ulimit -s)
    auto sys_ptr = std::make_unique<System<L>>();
    System<L>& sys = *sys_ptr;
    Measurements meas;

    std::vector<System<L>> systems;
    std::vector<Measurements> measurements;
    std::vector<int> map;  // key: replica ID, value: temperature ID
    std::vector<int> tmap; // key: temperature ID, value: replica ID

    std::vector<RealType> Results;

    for (int i = 0; i < 4; ++i) {
        sys.rng.state[i] = master_state[i];
    }

    meas.init_meas(mcs_in, pt_mcs_in, burnin_in, lag_in, Ti_in, Tf_in, n_rep);

    // some diagnostics to keep track of execution in a cluster
    std::cout << "L=" << L << "  temperatures: " << meas.n_rep << " on [" << Ti_in << ", " << Tf_in << "]";
    if (meas.n_pad > 0)
        std::cout << " + " << meas.n_pad << " above Tf (hottest T=" << 1.0 / meas.beta_of(0) << ")";
    std::cout << " = " << meas.length << " replicas on " << meas.Nthreads << " threads" << std::endl;
    sys.init_sys(random_start, J1, J2, J3, Jp);
    sys.set_neighbors();

    /* ONE disorder realization per call to Run(), drawn here and then copied into every replica below,
     * so all temperatures share the same sample. Different realizations come from different values of
     * int_param (different long_jump()s).
     */
    sys.set_nematic_disorder(J1, W, kernel);

    std::vector<RealType> slot_beta(std::size_t(meas.length));
    for (int k = 0; k < meas.length; ++k)
        slot_beta[std::size_t(k)] = meas.beta_of(k);
    System<L>::build_exp_tables(slot_beta);

    // Push the master stream past the block just consumed by the disorder, so the thermal streams created
    // below cannot overlap with it. 
    if (W != 0.0) {
        sys.rng.jump();
        for (int i = 0; i < 4; ++i)
            master_state[i] = sys.rng.state[i];
    }

    // For checking neighbor vector manually
    // for (size_t i = 0; i < sys.NbSize; ++i) {
    //     if (fast_mod(i,sys.Z) == 0) { // new site 
    //         std::cout << "---------------------------------------------" << std::endl;
    //         std::cout << "Neighbors of site = " << i/sys.Z << " are:\n";
    //     }
    //     std::cout << sys.Neighbors[i] << " ";
    // }

    systems.assign(meas.length, sys);
    measurements.assign(meas.length, meas);
    map.assign(meas.length, 0);
    tmap.assign(meas.length, 0);

    // 0:m2, 1:m4, 2:E, 3:E^2, 4:O^2, 5:O^4, 6:Nl^2, 7:Nl^4, 8:Ng^2, 9:Ng^4, 10:|psi_1|, 11:|psi_2|,
    // 12:|psi_3|, 13:|psi_1 psi_2 psi_3|, 14:eta_1, 15:eta_2, 16:eta_3, 17:|psi|, 18:|Nl|, 19:|Ng|.
    Results.assign(20 * meas.length, 0.0);

    /* parallel-tempering diagnostics 
     * swap_try/swap_ok[k] count attempts and acceptances for the pair of temperature slots (k, k+1).
     * flow[replica] records which end of the chain that replica visited last (+1 cold, -1 hot); a
     * round trip is counted whenever a replica reaches the cold end having previously reached the
     * hot end. 
     *
     * Note that high pairwise acceptance is necessary but not sufficient, since a single bottleneck
     * stops replicas traversing the chain...
     */
    std::vector<IntType> swap_try(meas.length - 1, 0), swap_ok(meas.length - 1, 0);
    std::vector<int8_t>  flow(meas.length, 0);
    IntType round_trips = 0;

    /* total-energy histograms 
     * One histogram of the TOTAL energy per temperature slot. 
     *
     * The energy is continuous once the disorder is switched on, so the bin range is measured in
     * a single pass just after burn-in and then padded. Samples that still fall outside are
     * counted separately rather than piled into the edge bins, where they would look like real
     * peaks. Set n_ebins = 0 in the input file to switch the whole thing off.
     */
    const int nbins = (n_ebins > 0) ? n_ebins : 0;
    std::vector<IntType> ehist, ehist_under, ehist_over;
    if (nbins > 0) {
        ehist.assign(std::size_t(meas.length) * nbins, 0);
        ehist_under.assign(meas.length, 0);
        ehist_over.assign(meas.length, 0);
    }
    RealType e_lo = 0.0, e_hi = 0.0, inv_bin = 0.0;
    RealType e_scan_lo =  std::numeric_limits<RealType>::max();
    RealType e_scan_hi = -std::numeric_limits<RealType>::max();

    /* spin configuration snapshots 
     * Every snap_every parallel-tempering rounds, the configuration of whichever replica is sitting
     * at the COLDEST temperature at that moment is appended to <o>_snap.dat. Snapshots are taken in
     * the measurement loop only, so the burn-in is already behind them; the same PT diagnostics that
     * validate the averages (round trips per replica in <o>_pt.dat) apply to them as well.
     *
     * (snap_every = 0 switches the snapshots off.)
     */
    const int snap_period = (snap_every > 0) ? snap_every : 0;
    std::string snap_file = output_file_name;
    if (snap_period > 0) {
        if (snap_file.size() > 4 && snap_file.compare(snap_file.size() - 4, 4, ".dat") == 0)
            snap_file = snap_file.substr(0, snap_file.size() - 4);
        snap_file += "_snap.dat";

        const int rounds_est = (mcs_in + (pt_mcs_in - mcs_in % pt_mcs_in) % pt_mcs_in) / pt_mcs_in;
        const int n_snaps    = rounds_est / snap_period;

        std::cout << "L=" << L << "  snapshots: " << n_snaps << " configurations to " << snap_file
                  << " (one every " << snap_period * pt_mcs_in * lag_in << " MCS, about " << std::endl;

        if (n_snaps == 0)
            std::cout << "   -> none will be written: snap_every=" << snap_period
                      << " exceeds the " << rounds_est << " swap rounds in this run." << std::endl;
    }

    // make sure each replica generates different random numbers
    for (int i = 0; i < meas.length; ++i) {
        for (int temp = 0; temp < 4; ++temp)
           systems[i].rng.state[temp] = master_state[temp];

        systems[i].rng.jump();
        // master state is always the deepest state in the prng sequence
        for (int j = 0; j < 4; ++j)
            master_state[j] = systems[i].rng.state[j];

        map[i]  = i;
        tmap[i] = i;
        measurements[i].current_beta = measurements[i].beta_of(i);
    }

    // A hot start must be drawn per replica, otherwise all replicas share one configuration
    if (random_start) {
        for (int i = 0; i < meas.length; ++i)
            for (int site = 0; site < System<L>::Total; ++site)
                systems[i].put(site, systems[i].rng.next_double() < 0.5);
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        // Equilibrate ALL replicas with parallel tempering
        int burnin_total = burnin_in + (pt_mcs_in - burnin_in%pt_mcs_in)%pt_mcs_in;
        for (int i = 0; i < (burnin_total / pt_mcs_in); ++i) {
            for (int pt_step = 0; pt_step < pt_mcs_in; ++pt_step) {
                for (int j = 0; j < (meas.length / meas.Nthreads); ++j) {
                    int idx = tid + j * meas.Nthreads;
                    for (int k = 0; k < lag_in; ++k) {
                        systems[idx].MCS_Metropolis(map[idx]); // one MCS
                    }
                }
            }

            for (int j = 0; j < (meas.length / meas.Nthreads); ++j) {
                int idx = tid + j * meas.Nthreads;
                Get(systems[idx], measurements[idx]); // collect observables
            }

            // start proposing swaps
            #pragma omp barrier 

            #pragma omp for schedule(static,1)
            for (int k = 0; k < meas.length - 1; k += 2) { // propose swaps for all even replicas
                RealType Delta_E = (measurements[tmap[k+1]].E1 - measurements[tmap[k]].E1) * System<L>::Total;
                RealType Delta_B = measurements[tmap[k+1]].current_beta - measurements[tmap[k]].current_beta;
                RealType exponent = Delta_B * Delta_E;

                if (exponent >= 0.0 || systems[tmap[k]].rng.next_double() < std::exp(exponent)) { // swap
                    RealType temp1                       = measurements[tmap[k]].current_beta;
                    measurements[tmap[k]].current_beta   = measurements[tmap[k+1]].current_beta;
                    measurements[tmap[k+1]].current_beta = temp1;
                    int temp2      = tmap[k];
                    tmap[k]        = tmap[k+1];
                    tmap[k+1]      = temp2;
                    map[tmap[k]]   = k;
                    map[tmap[k+1]] = k+1;
                }
            }

            #pragma omp for schedule(static,1)
            for (int k = 1; k < meas.length - 1; k += 2) { // propose swaps for all odd replicas
                RealType Delta_E = (measurements[tmap[k+1]].E1 - measurements[tmap[k]].E1) * System<L>::Total;
                RealType Delta_B = measurements[tmap[k+1]].current_beta - measurements[tmap[k]].current_beta;
                RealType exponent = Delta_B * Delta_E;

                if (exponent >= 0.0 || systems[tmap[k]].rng.next_double() < std::exp(exponent)) { // swap
                    RealType temp1                       = measurements[tmap[k]].current_beta;
                    measurements[tmap[k]].current_beta   = measurements[tmap[k+1]].current_beta;
                    measurements[tmap[k+1]].current_beta = temp1;
                    int temp2      = tmap[k];
                    tmap[k]        = tmap[k+1];
                    tmap[k+1]      = temp2;
                    map[tmap[k]]   = k;
                    map[tmap[k+1]] = k+1;
                }
            }
        } // end of burn in phase

        /* One pass over every replica to find the span of the total energy, so the histogram
         * bins cover what is actually sampled instead of the (far wider) theoretical range.
         */
        if (nbins > 0) {
            #pragma omp for schedule(static,1) reduction(min:e_scan_lo) reduction(max:e_scan_hi)
            for (int i = 0; i < meas.length; ++i) {
                Get(systems[i], measurements[i]);
                RealType e = measurements[i].E1 * RealType(System<L>::Total);
                if (e < e_scan_lo)
                    e_scan_lo = e;
                if (e > e_scan_hi)
                    e_scan_hi = e;
            }
            #pragma omp single
            {
                RealType span = e_scan_hi - e_scan_lo;
                if (!(span > 0.0))
                    span = 1.0;
                e_lo    = e_scan_lo - 0.20 * span;   // room for drift and fluctuation beyond burn-in
                e_hi    = e_scan_hi + 0.20 * span;
                inv_bin = RealType(nbins) / (e_hi - e_lo);
            }
        }

        // make total_mcs a multiple of pt_mcs
        int total_mcs = mcs_in + (pt_mcs_in - mcs_in%pt_mcs_in)%pt_mcs_in;
        for (int i = 0; i < (total_mcs / pt_mcs_in); ++i) {
            for (int pt_step = 0; pt_step < pt_mcs_in; ++pt_step) {
                /* Evolve and measure every replica this thread owns.
                 *
                 * Each observable is accumulated straight into the slot of the temperature the
                 * replica currently sits at. map is a bijection replica -> temperature, and the
                 * replicas are partitioned across threads by idx, so no two threads ever touch the
                 * same entry of Results (or of the histogram) and no synchronization is needed here.
                 */
                for (int j = 0; j < (meas.length / meas.Nthreads); ++j){
                    int idx = tid + j * meas.Nthreads;

                    for (int l = 0; l < lag_in; ++l) {
                        systems[idx].MCS_Metropolis(map[idx]);
                    }
                    Get(systems[idx], measurements[idx]);

                    const int k_temp = map[idx];
                    const Measurements& mi = measurements[idx];

                    Results[0*meas.length + k_temp] += mi.m2_cs;
                    Results[1*meas.length + k_temp] += mi.m4_cs;
                    Results[2*meas.length + k_temp] += mi.E1;
                    Results[3*meas.length + k_temp] += mi.E1 * mi.E1;
                    Results[4*meas.length + k_temp] += mi.O2_cs;
                    Results[5*meas.length + k_temp] += mi.O4_cs;
                    Results[6*meas.length + k_temp] += mi.Nl2_cs;
                    Results[7*meas.length + k_temp] += mi.Nl4_cs;
                    Results[8*meas.length + k_temp] += mi.Ng2_cs;
                    Results[9*meas.length + k_temp] += mi.Ng4_cs;

                    Results[10*meas.length + k_temp] += mi.psi_abs_cs[0];
                    Results[11*meas.length + k_temp] += mi.psi_abs_cs[1];
                    Results[12*meas.length + k_temp] += mi.psi_abs_cs[2];

                    Results[13*meas.length + k_temp] += mi.psi_prod_cs;

                    Results[14*meas.length + k_temp] += mi.eta_cs[0];
                    Results[15*meas.length + k_temp] += mi.eta_cs[1];
                    Results[16*meas.length + k_temp] += mi.eta_cs[2];

                    Results[17*meas.length + k_temp] += mi.psi_norm_cs;

                    Results[18*meas.length + k_temp] += mi.Nl_abs_cs;
                    Results[19*meas.length + k_temp] += mi.Ng_abs_cs;

                    if (nbins > 0) {
                        const RealType e = mi.E1 * RealType(System<L>::Total);
                        int b = int(std::floor((e - e_lo) * inv_bin));
                        if      (b < 0)      ++ehist_under[k_temp];
                        else if (b >= nbins) ++ehist_over[k_temp];
                        else                 ++ehist[std::size_t(k_temp) * nbins + b];
                    }
                }
            }

            // Every replica's energy must be in place before the swaps read it
            #pragma omp barrier

            // Parallel tempering swaps
            #pragma omp for schedule(static,1)
            for (int k = 0; k < meas.length - 1; k += 2) {
                // Multiply by Total to convert energy-per-spin back to raw energy for swap acceptance
                RealType Delta_E = (measurements[tmap[k+1]].E1 - measurements[tmap[k]].E1) * System<L>::Total;
                RealType Delta_B = measurements[tmap[k+1]].current_beta - measurements[tmap[k]].current_beta;
                RealType exponent = Delta_B * Delta_E;

                ++swap_try[k];

                if (exponent >= 0.0 || systems[tmap[k]].rng.next_double() < std::exp(exponent)) {
                    ++swap_ok[k];
                    RealType temp1                       = measurements[tmap[k]].current_beta;
                    measurements[tmap[k]].current_beta   = measurements[tmap[k+1]].current_beta;
                    measurements[tmap[k+1]].current_beta = temp1;
                    int temp2      = tmap[k];
                    tmap[k]        = tmap[k+1];
                    tmap[k+1]      = temp2;
                    map[tmap[k]]   = k;
                    map[tmap[k+1]] = k+1;
                }
            }

            #pragma omp for schedule(static,1)
            for (int k = 1; k < meas.length - 1; k += 2) {
                RealType Delta_E = (measurements[tmap[k+1]].E1 - measurements[tmap[k]].E1) * System<L>::Total;
                RealType Delta_B = measurements[tmap[k+1]].current_beta - measurements[tmap[k]].current_beta;
                RealType exponent = Delta_B * Delta_E;

                ++swap_try[k];

                if (exponent >= 0.0 || systems[tmap[k]].rng.next_double() < std::exp(exponent)) {
                    ++swap_ok[k];
                    RealType temp1                       = measurements[tmap[k]].current_beta;
                    measurements[tmap[k]].current_beta   = measurements[tmap[k+1]].current_beta;
                    measurements[tmap[k+1]].current_beta = temp1;
                    int temp2      = tmap[k];
                    tmap[k]        = tmap[k+1];
                    tmap[k+1]      = temp2;
                    map[tmap[k]]   = k;
                    map[tmap[k+1]] = k+1;
                }
            }

            /* Round trips. The omp for above ends with an implicit barrier, so tmap is settled here.
             * A replica that reaches the cold end after having reached the hot end has completed one
             * traversal of the whole temperature chain.
             */
            #pragma omp single
            {
                /* Temperature slot k carries beta = Bf + dB*k, so slot 0 is the HIGH temperature end
                 * and slot length-1 the LOW one: tmap[0] is the replica currently at the hot end and
                 * tmap[length-1] the one at the cold end. 
                 */
                // (with padding above Tf, slot k carries beta_of(k); slot 0 is still the hot end)
                const int hot  = tmap[0];
                const int cold = tmap[meas.length - 1];
                if (flow[cold] == -1)
                    ++round_trips;
                flow[cold] = +1;
                flow[hot]  = -1;

                /* Configuration snapshot of whichever replica is at the lowest temperature right now.
                 */
                if (snap_period > 0 && (i + 1) % snap_period == 0) {
                    const Measurements& mc = measurements[cold];

                    std::ofstream snap(snap_file, std::ios::app);
                    snap << "# snapshot L=" << L
                         << " round=" << (i + 1)
                         << " mcs=" << IntType(i + 1) * pt_mcs_in * lag_in
                         << " T=" << std::setprecision(10) << 1.0 / mc.current_beta
                         << " replica=" << cold
                         << " E=" << mc.E1
                         << " O2=" << mc.O2_cs
                         << " psi_abs=" << mc.psi_abs_cs[0] << ","
                                        << mc.psi_abs_cs[1] << ","
                                        << mc.psi_abs_cs[2]
                         << std::setprecision(6) << "\n";
                    snap << "# " << L * L << " lines of " << 2 * L << " characters: line k*L+i,"
                            " column 2*j+s is the spin of cell (i,j) in layer k on sublattice s"
                            " ('+' up, '-' down), i.e. site index 2*((k*L+i)*L+j)+s\n";

                    /* One line per (layer, row). Those 2L sites are contiguous in the site index,
                     * so the line is just a slice of Status.
                     */
                    const System<L>& sc = systems[cold];
                    std::string row;
                    row.reserve(std::size_t(2 * L));
                    for (int site = 0; site < System<L>::Total; ++site) {
                        row.push_back(sc.spin(site) > 0 ? '+' : '-');
                        if (row.size() == std::size_t(2 * L)) {
                            snap << row << "\n";
                            row.clear();
                        }
                    }
                    snap << "\n";
                }
            }
        }

        #pragma omp single
        {
            RealType norm = 1.0 / RealType(total_mcs);
            // RealType U_O, U_sg, Cv;

            /* Only the n_rep temperatures on [Ti, Tf] are written, so this file has the same rows
             * for every L and whatever the thread count. The padding slots above Tf are still simulated
             * and still appear in <o>_pt.dat and <o>_ehist.dat.
             * (To write them here too, start at k = 0.)
             */
            for (int k = meas.n_pad; k < meas.length; ++k) {
                RealType m2 = Results[0*meas.length + k] * norm;
                RealType m4 = Results[1*meas.length + k] * norm;
                RealType e1 = Results[2*meas.length + k] * norm;
                RealType e2 = Results[3*meas.length + k] * norm;
                RealType O2 = Results[4*meas.length + k] * norm;
                RealType O4 = Results[5*meas.length + k] * norm;

                RealType Nl2 = Results[6*meas.length + k] * norm;
                RealType Nl4 = Results[7*meas.length + k] * norm;
                RealType Ng2 = Results[8*meas.length + k] * norm;
                RealType Ng4 = Results[9*meas.length + k] * norm;

                RealType psi1_abs = Results[10*meas.length + k] * norm;
                RealType psi2_abs = Results[11*meas.length + k] * norm;
                RealType psi3_abs = Results[12*meas.length + k] * norm;

                RealType psi_prod = Results[13*meas.length + k] * norm;

                RealType eta1_av = Results[14*meas.length + k] * norm;
                RealType eta2_av = Results[15*meas.length + k] * norm;
                RealType eta3_av = Results[16*meas.length + k] * norm;

                RealType psi_norm = Results[17*meas.length + k] * norm;  

                RealType Nl_abs = Results[18*meas.length + k] * norm;   
                RealType Ng_abs = Results[19*meas.length + k] * norm;   

                RealType beta = measurements[tmap[k]].current_beta;

                std::ofstream out_append(output_file_name, std::ios::app);
                out_append << std::scientific << std::setprecision(12) << L << " "
                        << 1.0 / beta << " "
                        << m2 << " " << m4 << " "
                        << e1 << " " << e2 << " "
                        << O2 << " " << O4 << " "
                        << Nl2 << " " << Nl4 << " "
                        << Ng2 << " " << Ng4 << " "
                        << psi1_abs << " " << psi2_abs << " " << psi3_abs << " "
                        << psi_prod << " "
                        << eta1_av << " " << eta2_av << " " << eta3_av << " "
                        << psi_norm << " "
                        << Nl_abs << " " << Ng_abs << "\n";
            }

            /* parallel-tempering report
             * Written to <output>_pt.dat, one block per L, and summarized on stdout.
             */
            std::string pt_file = output_file_name;
            if (pt_file.size() > 4 && pt_file.compare(pt_file.size() - 4, 4, ".dat") == 0)
                pt_file = pt_file.substr(0, pt_file.size() - 4);
            pt_file += "_pt.dat";

            std::ofstream pt_out(pt_file, std::ios::app);
            pt_out << "# L=" << L << " temperatures=" << meas.length
                   << " n_rep=" << meas.n_rep << " padding=" << meas.n_pad;
            if (meas.n_pad > 0)
                pt_out << " (slots 0.." << meas.n_pad - 1 << " are above Tf)";
            pt_out << " pt_mcs=" << pt_mcs_in << "\n";
            pt_out << "# k / T_k / T_k+1 / attempts / accepted / rate / dbeta*sigma_E\n";

            RealType rate_min = 1.0, rate_sum = 0.0, x_max = 0.0;
            int      k_worst  = 0;

            for (int k = 0; k < meas.length - 1; ++k) {
                RealType T_k  = 1.0 / meas.beta_of(k);
                RealType T_k1 = 1.0 / meas.beta_of(k + 1);
                RealType rate = swap_try[k] > 0 ? RealType(swap_ok[k]) / RealType(swap_try[k]) : 0.0;

                /* Independent cross-check on the measured rate: the swap acceptance is controlled by
                 * dbeta * sigma_E, with sigma_E = N * sqrt(<e^2> - <e>^2) the width of the total-energy
                 * distribution. Values around 1.5 give a usable acceptance; much above that and the
                 * chain is severed here.
                 */
                RealType e1  = Results[2 * meas.length + k] * norm;
                RealType e2  = Results[3 * meas.length + k] * norm;
                RealType var = e2 - e1 * e1;
                RealType x   = std::fabs(meas.beta_of(k + 1) - meas.beta_of(k)) * RealType(System<L>::Total)
                             * std::sqrt(var > 0.0 ? var : 0.0);

                if (rate < rate_min) { rate_min = rate; k_worst = k; }
                rate_sum += rate;
                if (x > x_max) x_max = x;

                pt_out << k << " " << T_k << " " << T_k1 << " "
                       << swap_try[k] << " " << swap_ok[k] << " "
                       << rate << " " << x << "\n";
            }

            RealType rate_mean = rate_sum / RealType(meas.length - 1);
            const IntType pt_rounds = total_mcs / pt_mcs_in;
            RealType trips_per_replica = RealType(round_trips) / RealType(meas.length);

            const char* verdict = "OK";
            if (rate_min < 0.05)
                verdict = "BROKEN (grid too coarse)";
            else if (trips_per_replica < 1.0)
                verdict = "SHORT (chain intact, run too short)";

            pt_out << "# summary L=" << L
                   << " rate_min=" << rate_min << " (pair " << k_worst << ")"
                   << " rate_mean=" << rate_mean
                   << " pt_rounds=" << pt_rounds
                   << " round_trips=" << round_trips
                   << " per_replica=" << trips_per_replica
                   << " max_dbeta_sigmaE=" << x_max
                   << " verdict=" << verdict << "\n\n";

            std::cout << "PT check L=" << L
                      << ": acceptance min=" << std::fixed << std::setprecision(3) << rate_min
                      << " mean=" << rate_mean
                      << ", round trips/replica=" << std::setprecision(2) << trips_per_replica
                      << " over " << pt_rounds << " swap rounds"
                      << ", " << verdict << std::defaultfloat << std::endl;

            /* total-energy histograms 
             * Written to <output>_ehist.dat, one block per L. Row k holds the histogram of the
             * total energy at temperature slot k. Bin b covers
             *     E_lo + b*(E_hi-E_lo)/nbins  ..  E_lo + (b+1)*(E_hi-E_lo)/nbins
             */
            if (nbins > 0) {
                std::string eh_file = output_file_name;
                if (eh_file.size() > 4 && eh_file.compare(eh_file.size() - 4, 4, ".dat") == 0)
                    eh_file = eh_file.substr(0, eh_file.size() - 4);
                eh_file += "_ehist.dat";

                IntType n_under = 0, n_over = 0;
                for (int k = 0; k < meas.length; ++k) { n_under += ehist_under[k]; n_over += ehist_over[k]; }
                const IntType n_total = IntType(total_mcs) * meas.length;

                std::ofstream eh(eh_file, std::ios::app);
                eh << "# ehist L=" << L << " temperatures=" << meas.length
                   << " padding=" << meas.n_pad
                   << " nbins=" << nbins
                   << " E_lo=" << std::setprecision(10) << e_lo
                   << " E_hi=" << e_hi
                   << " samples_per_T=" << total_mcs
                   << " outside=" << (n_under + n_over)
                   << " outside_frac=" << RealType(n_under + n_over) / RealType(n_total)
                   << std::setprecision(6) << "\n";
                eh << "# k / T / underflow / overflow / counts[0..nbins-1]\n";

                for (int k = 0; k < meas.length; ++k) {
                    eh << k << " " << 1.0 / meas.beta_of(k) << " "
                       << ehist_under[k] << " " << ehist_over[k];
                    for (int b = 0; b < nbins; ++b)
                        eh << " " << ehist[std::size_t(k) * nbins + b];
                    eh << "\n";
                }
                eh << "\n";
            }
        } // omp single
    } // omp parallel
} // Run()


// Compiler routes runtime parameters to target template instantiations
void dispatch_L(int L, int mcs, int pt_mcs, int burnin, int lag,
                RealType Ti, RealType Tf, RealType J1, RealType J2, RealType J3, RealType Jp,
                RealType W, int kernel, int n_rep, int n_ebins, int snap_every,
                bool random_start, const std::string& output) {
    switch(L) {
        case 6:
            Run<6>(mcs, pt_mcs, burnin, lag, Ti, Tf, J1, J2, J3, Jp, W, kernel, n_rep, n_ebins, snap_every, random_start, output);
            break;
        case 12:
            Run<12>(mcs, pt_mcs, burnin, lag, Ti, Tf, J1, J2, J3, Jp, W, kernel, n_rep, n_ebins, snap_every, random_start, output);
            break;
        case 24:
            Run<24>(mcs, pt_mcs, burnin, lag, Ti, Tf, J1, J2, J3, Jp, W, kernel, n_rep, n_ebins, snap_every, random_start, output);
            break;
        case 48:
            Run<48>(mcs, pt_mcs, burnin, lag, Ti, Tf, J1, J2, J3, Jp, W, kernel, n_rep, n_ebins, snap_every, random_start, output);
            break;
        default:
            std::cerr << "Error: System size L=" << L << " is not supported as a template parameter." << std::endl;
            break;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Please provide a positive integer to seed the PRNG." << std::endl;
        return 1;
    }

    int int_param = std::stoi(argv[1]);
    int init_L, n_sizes, mcs, pt_mcs, burnin, lag, seed, n_ebins, kernel, snap_every, n_rep;
    bool random_start;
    RealType Ti, Tf, J1, J2, J3, Jp, W;
    Xoshiro256Plus initial_prng;

    /* NOTE:
     *
     *   W        : nematic random-field strength (W = 0 recovers the clean model exactly)
     *   n_ebins  : bins in the per-temperature total-energy histogram written to <output>_ehist.dat.
     *              (0 switches the histograms off.)
     *   kernel   : spatial correlation kernel for the nematic random field.
     *              0 = uncorrelated (identical to the previous behaviour), 1 = 6-site, 2 = 10-site,
     *              3 = 13-site, 4 = 18-site, 5 = 35-site. 
     *   snap_every: how often to append a configuration of the coldest replica to <o>_snap.dat,
     *              counted in parallel-tempering rounds.
     *              (0 switches them off.)
     *   n_rep     : number of replicas in [Ti, Tf]. Same for every system size L.
     *
     * Input file layout (all 19 values required, in this order):
     *
     *   init_L n_sizes random_start Ti Tf J1 J2 J3 Jp mcs pt_mcs burnin lag seed W n_ebins kernel snap_every n_rep
     *
     * N.B.: a decimal in an integer field (e.g. 2.5 for n_ebins) is an error rather than being cut to 2.
     */
    static const char* input_names[] = { "init_L", "n_sizes", "random_start", "Ti", "Tf", "J1", "J2", "J3", "Jp",
                                         "mcs", "pt_mcs", "burnin", "lag", "seed", "W", "n_ebins", "kernel",
                                         "snap_every", "n_rep" };
    constexpr std::size_t n_input = sizeof(input_names) / sizeof(input_names[0]);

    std::ifstream infile("NNIH_input.in");
    if (!infile) {
        std::cerr << "Error: cannot open NNIH_input.in" << std::endl;
        return 1;
    }

    std::vector<std::string> tokens;
    for (std::string t; infile >> t; )
        tokens.push_back(t);

    if (tokens.size() != n_input) {
        std::cerr << "Error: NNIH_input.in must hold exactly " << n_input << " values, found " << tokens.size()
                  << ". Expected order:\n  ";
        for (const char* name : input_names)
            std::cerr << " " << name;
        std::cerr << std::endl;
        return 1;
    }

    std::size_t next_token = 0;
    auto read_value = [&](auto& value) -> bool {
        std::istringstream is(tokens[next_token]);
        char leftover;
        const bool ok = static_cast<bool>(is >> value) && !(is >> leftover);
        if (!ok)
            std::cerr << "Error: cannot read " << input_names[next_token] << " from '" << tokens[next_token]
                      << "' (value " << next_token + 1 << " of " << n_input << ")." << std::endl;
        ++next_token;
        return ok;
    };

    if (!(read_value(init_L) && read_value(n_sizes) && read_value(random_start) && read_value(Ti) && read_value(Tf)
          && read_value(J1) && read_value(J2) && read_value(J3) && read_value(Jp) && read_value(mcs)
          && read_value(pt_mcs) && read_value(burnin) && read_value(lag) && read_value(seed) && read_value(W)
          && read_value(n_ebins) && read_value(kernel) && read_value(snap_every) && read_value(n_rep)))
        return 1;

    if (kernel < 0 || kernel > 5) {
        std::cerr << "Error: kernel must be 0, 1, 2, 3, 4 or 5 (got " << kernel << ")." << std::endl;
        return 1;
    }

    if (snap_every < 0) {
        std::cerr << "Error: snap_every must not be negative (got " << snap_every << ")." << std::endl;
        return 1;
    }

    if (n_rep < 2) {
        std::cerr << "Error: n_rep must be at least 2 (got " << n_rep << ")." << std::endl;
        return 1;
    }

    std::string output = "NNIH_" + std::to_string(int_param) + ".dat";

    // If this file exists, it's probably because that cluster node went offline and back on again
    // Its best not to overwrite or append anything to it to avoid corrupting the results...
    if (std::filesystem::exists(output)) {
        std::cout << "File: " << output << " already exists." << std::endl;
        return 1;
   }

    std::ofstream outfile(output);
    outfile << "# dim=3\n";
    outfile << "# J1=" << J1 << " J2=" << J2 << " J3=" << J3 << " Jp=" << Jp << "\n";
    outfile << "# W=" << W << " kernel=" << kernel << "\n";
    outfile << "# n_rep=" << n_rep << " n_ebins=" << n_ebins << " snap_every=" << snap_every << "\n";
    outfile << "# mcs=" << mcs << "\n";
    outfile << "# pt_mcs=" << pt_mcs << "\n";
    outfile << "# burnin=" << burnin << "\n";
    outfile << "# seed=" << seed << "\n";
    outfile << "# random start=" << random_start << "\n";
    outfile << "# L / T / <m^2> / <m^4> / <E> / <E^2> / <O^2> / <O^4>"
               " / <Nl^2> / <Nl^4> / <Ng^2> / <Ng^4>"
               " / <|psi_1|> / <|psi_2|> / <|psi_3|>"
               " / <|psi_1 psi_2 psi_3|>"
               " / <eta_1> / <eta_2> / <eta_3>"
               " / <|psi|> / <|Nl|> / <|Ng|> \n";
    outfile.close();

    // to find out the initial state of the prng at the seed, instantiate it and call seed()
    initial_prng.seed(seed);

   // Now set prng to a unique spot for this current node. (other nodes long_jump() accordingly until they reach a unique sequence)
    for (int i = 0; i < int_param; ++i) {
        initial_prng.long_jump();
    }

    // Now copy this initial state to the master state which will dictate from which part of the prng sequence this program will run on
    // (from now on, never call long_jump() again)
    for (int i = 0; i < 4; ++i) {
        master_state[i] = initial_prng.state[i];
    }

    int L = init_L;
    for (int irrelevant_var = 0; irrelevant_var < n_sizes; ++irrelevant_var) {
        dispatch_L(L, mcs, pt_mcs, burnin, lag, Ti, Tf, J1, J2, J3, Jp, W, kernel, n_rep, n_ebins, snap_every, random_start, output);
        L *= 2;
    }

    return 0;
}
