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
    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

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
/*
#include <armadillo>

float calcRho(float iHdist, float iVdist, float iLdist) const;
static std::string description(bool full=true);
std::string name() const {return "oi";};
bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
enum Type {TypeTemperature, TypePrecipitation};
enum TransformType {TransformTypeNone, TransformTypeBoxCox};
float mVLength;
float mHLength;
float mHLengthC;
float mMu;
float mGamma;
std::string mBiasVariable;
int mMaxLocations;
float mSigma;
float mSigmaC;
float mDelta;
std::string mDeltaVariable;
std::string mNumVariable;
float mC;
float mEpsilonC;
float mEpsilon;
int mX;
int mY;
float mMinRho;
float mMaxBytes;
int mMinValidEns;
bool mSaveDiff;
float mElevGradient;
bool mExtrapolate;
float mWMin;
float mNewDeltaVar;
float mMaxElevDiff;
bool mDiagnose;
bool mLandOnly;
std::string mDiaFile;
bool mUseEns;
typedef arma::mat mattype;
typedef arma::vec vectype;
typedef arma::cx_mat cxtype;
float mLambda;
bool mCrossValidate;
Type mType;
TransformType mTransformType;
float calcDelta(float iOldDelta, const vec2& iY) const;
float transform(float iValue) const;
float invTransform(float iValue) const;
*/
