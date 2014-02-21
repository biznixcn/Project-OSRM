#ifndef EDGE_BASED_NODE_H
#define EDGE_BASED_NODE_H

#include "../Util/MercatorUtil.h"
#include "../Util/SimpleLogger.h"
#include "../typedefs.h"

#include <osrm/Coordinate.h>

#include <boost/assert.hpp>

#include <limits>

struct EdgeBasedNode {

    EdgeBasedNode() :
        forward_edge_based_node_id(std::numeric_limits<int>::max()),
        reverse_edge_based_node_id(std::numeric_limits<int>::max()),
        lat1(std::numeric_limits<int>::max()),
        lon1(std::numeric_limits<int>::max()),
        lat2(std::numeric_limits<int>::max()),
        lon2(std::numeric_limits<int>::max() >> 1),
        belongsToTinyComponent(false),
        name_id(std::numeric_limits<unsigned>::max()),
        forward_weight(std::numeric_limits<int>::max() >> 1),
        reverse_weight(std::numeric_limits<int>::max() >> 1),
        forward_offset(0),
        reverse_offset(0)
    { }

    EdgeBasedNode(
        NodeID forward_edge_based_node_id,
        NodeID reverse_edge_based_node_id,
        int lat1,
        int lon1,
        int lat2,
        int lon2,
        bool belongsToTinyComponent,
        NodeID name_id,
        int forward_weight,
        int reverse_weight,
        int forward_offset,
        int reverse_offset
    ) :
        forward_edge_based_node_id(forward_edge_based_node_id),
        reverse_edge_based_node_id(reverse_edge_based_node_id),
        lat1(lat1),
        lon1(lon1),
        lat2(lat2),
        lon2(lon2),
        belongsToTinyComponent(belongsToTinyComponent),
        name_id(name_id),
        forward_weight(forward_weight),
        reverse_weight(reverse_weight),
        forward_offset(forward_offset),
        reverse_offset(reverse_offset)
    { }

    double ComputePerpendicularDistance(
        const FixedPointCoordinate& query_location,
        FixedPointCoordinate & nearest_location,
        double & r
    ) const {
        BOOST_ASSERT( query_location.isValid() );

        const double x = lat2y(query_location.lat/COORDINATE_PRECISION);
        const double y = query_location.lon/COORDINATE_PRECISION;
        const double a = lat2y(lat1/COORDINATE_PRECISION);
        const double b = lon1/COORDINATE_PRECISION;
        const double c = lat2y(lat2/COORDINATE_PRECISION);
        const double d = lon2/COORDINATE_PRECISION;
        double p,q/*,mX*/,nY;
        if( std::abs(a-c) > std::numeric_limits<double>::epsilon() ){
            const double m = (d-b)/(c-a); // slope
            // Projection of (x,y) on line joining (a,b) and (c,d)
            p = ((x + (m*y)) + (m*m*a - m*b))/(1. + m*m);
            q = b + m*(p - a);
        } else {
            p = c;
            q = y;
        }
        nY = (d*p - c*q)/(a*d - b*c);

        //discretize the result to coordinate precision. it's a hack!
        if( std::abs(nY) < (1./COORDINATE_PRECISION) ) {
            nY = 0.;
        }

        r = (p - nY*a)/c;// These values are actually n/m+n and m/m+n , we need
        // not calculate the explicit values of m an n as we
        // are just interested in the ratio
        if( std::isnan(r) ) {
            r = ((lat2 == query_location.lat) && (lon2 == query_location.lon)) ? 1. : 0.;
        } else if( std::abs(r) <= std::numeric_limits<double>::epsilon() ) {
            r = 0.;
        } else if( std::abs(r-1.) <= std::numeric_limits<double>::epsilon() ) {
            r = 1.;
        }
        BOOST_ASSERT( !std::isnan(r) );
        if( r <= 0. ){
            nearest_location.lat = lat1;
            nearest_location.lon = lon1;
        } else if( r >= 1. ){
            nearest_location.lat = lat2;
            nearest_location.lon = lon2;
        } else {
            // point lies in between
            nearest_location.lat = y2lat(p)*COORDINATE_PRECISION;
            nearest_location.lon = q*COORDINATE_PRECISION;
        }
        BOOST_ASSERT( nearest_location.isValid() );

        // TODO: Replace with euclidean approximation when k-NN search is done
        // const double approximated_distance = FixedPointCoordinate::ApproximateEuclideanDistance(
        const double approximated_distance = FixedPointCoordinate::ApproximateDistance(
            query_location,
            nearest_location
        );
        BOOST_ASSERT( 0. <= approximated_distance );
        return approximated_distance;
    }

    inline FixedPointCoordinate Centroid() const {
        FixedPointCoordinate centroid;
        //The coordinates of the midpoint are given by:
        //x = (x1 + x2) /2 and y = (y1 + y2) /2.
        centroid.lon = (std::min(lon1, lon2) + std::max(lon1, lon2))/2;
        centroid.lat = (std::min(lat1, lat2) + std::max(lat1, lat2))/2;
        return centroid;
    }

    bool IsCompressed() {
        return (forward_offset*reverse_offset) == 0;
    }

    NodeID forward_edge_based_node_id;
    NodeID reverse_edge_based_node_id;
    int lat1;
    int lon1;
    int lat2;
    int lon2:31;
    bool belongsToTinyComponent:1;
    NodeID name_id;
    int forward_weight;
    int reverse_weight;
    int forward_offset;
    int reverse_offset;
};

#endif //EDGE_BASED_NODE_H
