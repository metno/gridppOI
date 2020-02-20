#include "gridppOI.h"

gridppOI::KDTree::KDTree(const vec& lats, const vec& lons) {
    vec x, y, z;
    int nx = 3;
    int ny = 1;
    int normtype = 2;
    gridppOI::KDTree::convert_coordinates(lats, lons, x, y, z);

    alglib::real_2d_array a;
    a.setlength(lats.size(), 4);
    for(int i = 0; i < lats.size(); i++) {
        a[i][0] = x[i];
        a[i][1] = y[i];
        a[i][2] = z[i];
        a[i][3] = i;
    }
    alglib::kdtreebuild(a, nx, ny, normtype, mTree);
}

int gridppOI::KDTree::get_num_neighbours(float lat, float lon, float radius) {
    alglib::real_1d_array b = gridppOI::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius);
    return num;
}

ivec gridppOI::KDTree::get_neighbours(float lat, float lon, float radius) {
    alglib::real_1d_array b = gridppOI::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        ret[i] = ans[i][3];
    }
    return ret;

}

ivec gridppOI::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) {
    alglib::real_1d_array b = gridppOI::KDTree::ll2ar(lat, lon);
    int num = alglib::kdtreequeryrnn(mTree, b, radius);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num);
    for(int i = 0; i < num; i++) {
        ret[i] = ans[i][3];
    }

    alglib::real_1d_array rdist;
    alglib::kdtreequeryresultsdistances(mTree, rdist);
    distances.resize(num);
    for(int i = 0; i < num; i++) {
        distances[i] = rdist[i];
    }
    return ret;
}

ivec gridppOI::KDTree::get_closest_neighbours(float lat, float lon, int num) {

    alglib::real_1d_array b = gridppOI::KDTree::ll2ar(lat, lon);
    int num_found = alglib::kdtreequeryknn(mTree, b, num);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    ret.resize(num_found);
    for(int i = 0; i < num_found; i++) {
        ret[i] = ans[i][3];
    }

    return ret;
}
alglib::real_1d_array gridppOI::KDTree::ll2ar(float lat, float lon) {
    alglib::real_1d_array b;
    float x, y, z;
    gridppOI::KDTree::convert_coordinates(lat, lon, x, y, z);
    b.setlength(3);
    b[0] = x;
    b[1] = y;
    b[2] = z;
    return b;
}
int gridppOI::KDTree::get_nearest_neighbour(float lat, float lon) {
    alglib::real_1d_array b = gridppOI::KDTree::ll2ar(lat, lon);
    int num_found = alglib::kdtreequeryknn(mTree, b, 1);

    ivec ret;
    alglib::real_2d_array ans;
    alglib::kdtreequeryresultsxy(mTree, ans);
    return ans[0][3];
}
bool gridppOI::KDTree::convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords) {
    int N = lats.size();
    x_coords.resize(N);
    y_coords.resize(N);
    z_coords.resize(N);
    for(int i = 0; i < N; i++) {
        convert_coordinates(lats[i], lons[i], x_coords[i], y_coords[i], z_coords[i]);
    }

    return true;
}

bool gridppOI::KDTree::convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord) {

    float earth_radius = 6.37e6;
    double lonr = M_PI / 180 * lon;
    double latr = M_PI / 180 * lat;
    // std::cout << lon << " " << lat << std::endl;
    x_coord = std::cos(latr) * std::cos(lonr) * earth_radius;
    y_coord = std::cos(latr) * std::sin(lonr) * earth_radius;
    z_coord = std::sin(latr) * earth_radius;
    return true;
}
