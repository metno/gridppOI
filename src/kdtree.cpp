#include "gridppOI.h"

gridppOI::KDTree::KDTree(const vec& lats, const vec& lons) {
    mLats = lats;
    mLons = lons;
    vec x, y, z;
    gridppOI::KDTree::convert_coordinates(lats, lons, x, y, z);

    for(int i = 0; i < lats.size(); i++) {
        point p(x[i], y[i], z[i]);
        mTree.insert(std::make_pair(p, i));
    }
}

int gridppOI::KDTree::get_num_neighbours(float lat, float lon, float radius) {
    ivec indices = get_neighbours(lat, lon, radius);
    return indices.size();
}

ivec gridppOI::KDTree::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) {
    float x, y, z;
    gridppOI::KDTree::convert_coordinates(lat, lon, x, y, z);
    ivec indices = get_neighbours(lat, lon, radius);

    int num = indices.size();
    distances.resize(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        gridppOI::KDTree::convert_coordinates(mLats[indices[i]], mLons[indices[i]], x1, y1, z1);
        distances[i] = gridppOI::KDTree::calc_distance(x, y, z, x1, y1, z1);
    }

    return indices;
}

ivec gridppOI::KDTree::get_neighbours(float lat, float lon, float radius) {
    float x, y, z;
    gridppOI::KDTree::convert_coordinates(lat, lon, x, y, z);

    box bx(point(x - radius, y - radius, z - radius), point(x + radius, y + radius, z + radius));
    std::vector<value> results;
    mTree.query(boost::geometry::index::within(bx), std::back_inserter(results));
    int num = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num; i++) {
        float x1, y1, z1;
        gridppOI::KDTree::convert_coordinates(mLats[results[i].second], mLons[results[i].second], x1, y1, z1);
        float dist = gridppOI::KDTree::calc_distance(x, y, z, x1, y1, z1);
        if(dist > 0 && dist <= radius) {
            ret.push_back(results[i].second);
        }
    }
    return ret;
}

ivec gridppOI::KDTree::get_closest_neighbours(float lat, float lon, int num) {
    float x, y, z;
    gridppOI::KDTree::convert_coordinates(lat, lon, x, y, z);
    point p(x, y, z);

    struct is_not_equal {
        bool operator()(value const& v) const {
            float x0 = v.first.get<0>();
            float y0 = v.first.get<1>();
            float z0 = v.first.get<2>();
            return p.get<0>() != x0 || p.get<1>() != y0 || p.get<2>() != z0;
        };
        point p;
    };
    is_not_equal s;
    s.p = p;

    std::vector<value> results;
    mTree.query(boost::geometry::index::nearest(p, num) && boost::geometry::index::satisfies(s), std::back_inserter(results));
    int num_found = results.size();

    ivec ret;
    ret.reserve(num);
    for(int i = 0; i < num_found; i++) {
        float x1, y1, z1;
        gridppOI::KDTree::convert_coordinates(mLats[results[i].second], mLons[results[i].second], x1, y1, z1);
        float dist = gridppOI::KDTree::calc_distance(x, y, z, x1, y1, z1);
        std::cout << "Distance " << dist << std::endl;
        if(dist > 0) {
            ret.push_back(results[i].second);
        }
    }
    return ret;
}
int gridppOI::KDTree::get_nearest_neighbour(float lat, float lon) {
    return get_closest_neighbours(lat, lon, 1)[0];
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
float gridppOI::KDTree::calc_distance(float lat1, float lon1, float lat2, float lon2) {
    if(!(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360)) {
        // std::cout << " Cannot calculate distance, invalid lat/lon: (" << lat1 << "," << lon1 << ") (" << lat2 << "," << lon2 << ")";
        // std::cout << '\n';
    }
    if(lat1 == lat2 && lon1 == lon2)
        return 0;

    double lat1r = deg2rad(lat1);
    double lat2r = deg2rad(lat2);
    double lon1r = deg2rad(lon1);
    double lon2r = deg2rad(lon2);
    double radiusEarth = 6.378137e6;

    double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                   + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                   + sin(lat1r)*sin(lat2r);
    double dist = acos(ratio)*radiusEarth;
    return (float) dist;
}
float gridppOI::KDTree::calc_distance(float x0, float y0, float z0, float x1, float y1, float z1) {
    return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1) + (z0 - z1)*(z0 - z1));
}
float gridppOI::KDTree::deg2rad(float deg) {
   return (deg * M_PI / 180);
}
