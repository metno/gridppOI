#include <vector>
#include <set>
#include <string>
#include <armadillo>
#include <libalglib/interpolation.h>
typedef std::vector<std::vector<float> > vec2;
typedef std::vector<float> vec;
// typedef std::vector<float> fvec;
typedef std::vector<int> ivec;

namespace gridppOI {

    /** Optimal interpolation
     *  @param input 2D background field
     *  @param blats 2D matrix of latitudes [degrees]
     *  @param blons 2D matrix of longitudes [degrees]
     *  @param belevs 2D matrix of altitudes [m]
     *  @param blafs 2D matrix of land area fractions [1]
     *  @param pobs 1D vector of observations
     *  @param pci 1D vector of ci
     *  @param plats 1D vector of observation latitudes [degrees]
     *  @param plons 1D vector of observation longitudes [degrees]
     *  @param pelelvs 1D vector of observation altitudes [m]
     *  @param plafs 1D vector of observation land area fractions [1]
     *  @param minRho Localization cut-off rho value
     *  @param hlength Horizontal decorrelation length scale [m]
     *  @param vlength Vertical decorrelation length scale [m]
     *  @param wmin Minimum correlation for Land-sea
     *  @param maxElevDiff Maximum elevation difference between obs and background [m]
     *  @param landOnly Only include observations on land
     *  @param maxLocations Maximum obs to use in localization
     *  @param elevGradient Vertical temperature gradient used to move from grid to point [*C/m]
     *  @param epsilon Ratio of observation to background error
     *  @param output Write analysis field to this matrix
    **/
    bool optimal_interpolation(const vec2& input,
            const vec2& blats,
            const vec2& blons,
            const vec2& belevs,
            const vec2& blafs,
            const vec& pobs,  // gObs
            const vec& pci,   // gCi
            const vec& plats,
            const vec& plons,
            const vec& pelevs,
            const vec& plafs,
            float minRho,
            float hlength,
            float vlength,
            float wmin,
            float maxElevDiff,
            bool landOnly,
            int maxLocations,
            float elevGradient,
            float epsilon,
            vec2& output);

    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

    float getDistance(float lat1, float lon1, float lat2, float lon2, bool approx=false);
    bool isValid(float iValue);
    float deg2rad(float deg);
    float rad2deg(float rad);
    void debug(std::string string);
    void error(std::string string);
    float calcRho(float iHDist, float iVDist, float iLDist, float hlength, float vlength, float wmin);

    static float MV;
    static float pi;
    static double radiusEarth;

    void check_vec(vec2 input, int Y, int X);
    void check_vec(vec input, int S);

    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };

    class KDTree {
        public:
            KDTree(const vec& lats, const vec& lons);

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             * */
            int get_nearest_neighbour(float lat, float lon);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius);

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances);

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius);

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num);


            /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
             *  @param lats vector of latitudes [deg]
             *  @param lons vector of longitudes [deg]
             *  @param x_coords vector of x-coordinates [m]
             *  @param y_coords vector of y-coordinates [m]
             *  @param z_coords vector of z-coordinates [m]
             * */
            static bool convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords);

            /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
             *  @param lat latitude [deg]
             *  @param lon longitude [deg]
             *  @param x_coord x-coordinate [m]
             *  @param y_coord y-coordinate [m]
             *  @param z_coord z-coordinate [m]
             * */
            static bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord);

        private:
            alglib::kdtree mTree;
            static alglib::real_1d_array ll2ar(float lat, float lon);

    };
}
