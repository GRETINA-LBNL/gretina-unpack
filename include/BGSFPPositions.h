#ifndef __BGSFPPOSITIONS_H
#define __BGSFPPOSITIONS_H

#include "Structures.h"

/* Define the geometry and strip numbers for the FP
   (x,y,z) positions of some corners of the active area for
   top(t), west(w) and east(e) DSSD chips are needed. 
   Outer point is furthest from the apex of the vacuum window, 
   early point is counter-clockwise from center of chip,
   and later point is clockwise from center of chip.
   The origin has been defined as the outside apex of the pyramid-like
   vacuum window.  The coordinate system orientation is as
   expected when looking at the back of the detector.  
   
   Positive x is the to the right, positive y is up and positive z
   is coming at you. 

   Looking at the detector from the back:
                Top strip 0  /\  Back strip 0
        Top front strips    /  \     Top back strips
	 go this way \\    /    \    go this way //
            Top strip 31  /   T  \  Back strip 31
	   Back strip 31 |\      /| West strip 31
     East back strips    | \    / |   West front strips
      go this way \\     |  \  /  |    go this way //
           Back strip 0  | E \/ W | West strip 0
           East strip 0 \    |   / Back strip 0
     East front strips    \   |  /    West back strips
      go this way ||       \  | /      go this way ||
             East strip 31  \ |/  Back strip 31
	      Group1 strip 0 \/  Group1 strip 31
                Group1 strips go this way ||
                                                                          */
struct rect_point xyz(int chip, int sf, int sb, int ran);

struct rect_point t_outer = {-0.282843,  5.290898, -9.422356}; /* use as start 
								  of top front,
								  top back */
struct rect_point t_later = { 4.242641,  2.678109, -5.727315}; /* use as end of
								  top back */
struct rect_point t_early = {-4.808326,  2.678109, -5.727315}; /* use as end 
								  of top 
								  front */


struct rect_point w_outer = { 4.723473, -2.400500, -9.422356}; /* use as start
								  of west
								  front, west
								  back */
struct rect_point w_later = { 0.197990, -5.013289, -5.727315}; /* use as end of
								  west back */
struct rect_point w_early = { 4.723473,  2.825078, -5.727315}; /* use as end of
								  west front */


struct rect_point e_outer = {-4.440631, -2.890398, -9.422356}; /* use as start
								  of east
								  front, east
								  back */
struct rect_point e_later = {-4.440631,  2.335180, -5.727315}; /* use as end of
								  east back */
struct rect_point e_early = { 0.084853, -5.503187, -5.727315}; /* use as end of
								  east front */


/* Now calculated the change in x, y, and z with increasing strip number. */
//struct rect_point tf_delta = { (t_early.x - w_outer.x)/32.,
//			       (t_early.y - w_outer.y)/32.,
//			       (t_early.z - w_outer.z)/32. }; /* Top front */
//struct rect_point tb_delta = { (t_later.x - w_outer.x)/32.,
//			       (t_later.y - w_outer.y)/32.,
//			       (t_later.z - w_outer.z)/32. }; /* Top back */

//struct rect_point wf_delta = { (w_early.x - w_outer.x)/32.,
//			       (w_early.y - w_outer.y)/32.,
//			       (w_early.z - w_outer.z)/32. }; /* West front */
//struct rect_point wb_delta = { (w_later.x - w_outer.x)/32.,
//			       (w_later.y - w_outer.y)/32.,
//			       (w_later.z - w_outer.z)/32. }; /* West back */

//struct rect_point ef_delta = { (e_early.x - e_outer.x)/32.,
//			       (e_early.y - e_outer.y)/32.,
//			       (e_early.z - e_outer.z)/32. }; /* Top front */
//struct rect_point eb_delta = { (e_later.x - e_outer.x)/32.,
//			       (e_later.y - e_outer.y)/32.,
//			       (e_later.z - e_outer.z)/32. }; /* Top back */

/* To get to the center of the pixel, the origin should be offset from the
   edge of the active area to the center of the appropriate corner pixel,
   that is half a pixel in each direction. */
//struct rect_point t_origin = { (t_outer.x + tf_delta.x/2. + tb_delta.x/2.),
//			       (t_outer.y + tf_delta.y/2. + tb_delta.y/2.),
//			       (t_outer.z + tf_delta.z/2. + tb_delta.z/2.) };
//struct rect_point w_origin = { (w_outer.x + wf_delta.x/2. + wb_delta.x/2.),
//			       (w_outer.y + wf_delta.y/2. + wb_delta.y/2.),
//			       (w_outer.z + wf_delta.z/2. + wb_delta.z/2.) };
//struct rect_point e_origin = { (e_outer.x + ef_delta.x/2. + eb_delta.x/2.),
//			       (e_outer.y + ef_delta.y/2. + eb_delta.y/2.),
//			       (e_outer.z + ef_delta.z/2. + eb_delta.z/2.) };

/* New calculations, because Ken thinks Mesytec renumbered the strips from 
   the front in reverse order.  So he redefines the origins and deltas for
   detectors. */
struct rect_point tf_delta = { (t_early.x - t_outer.x)/32.,
			       (t_early.y - t_outer.y)/32.,
			       (t_early.z - t_outer.z)/32. }; /* Top front */
struct rect_point tb_delta = { (t_later.x - t_outer.x)/32.,
			       (t_later.y - t_outer.y)/32.,
			       (t_later.z - t_outer.z)/32. }; /* Top back */

struct rect_point wf_delta = { (w_early.x - w_outer.x)/32.,
			       (w_early.y - w_outer.y)/32.,
			       (w_early.z - w_outer.z)/32. }; /* West front */
struct rect_point wb_delta = { (w_later.x - w_outer.x)/32.,
			       (w_later.y - w_outer.y)/32.,
			       (w_later.z - w_outer.z)/32. }; /* West back */

struct rect_point ef_delta = { (e_early.x - e_outer.x)/32.,
			       (e_early.y - e_outer.y)/32.,
			       (e_early.z - e_outer.z)/32. }; /* Top front */
struct rect_point eb_delta = { (e_later.x - e_outer.x)/32.,
			       (e_later.y - e_outer.y)/32.,
			       (e_later.z - e_outer.z)/32. }; /* Top back */

/* To get to the center of the pixel, the origin should be offset from the
   edge of the active area to the center of the appropriate corner pixel,
   that is half a pixel in each direction. */
struct rect_point t_origin = { (t_outer.x + tf_delta.x/2. + tb_delta.x/2.),
			       (t_outer.y + tf_delta.y/2. + tb_delta.y/2.),
			       (t_outer.z + tf_delta.z/2. + tb_delta.z/2.) };
struct rect_point w_origin = { (w_outer.x + wf_delta.x/2. + wb_delta.x/2.),
			       (w_outer.y + wf_delta.y/2. + wb_delta.y/2.),
			       (w_outer.z + wf_delta.z/2. + wb_delta.z/2.) };
struct rect_point e_origin = { (e_outer.x + ef_delta.x/2. + eb_delta.x/2.),
			       (e_outer.y + ef_delta.y/2. + eb_delta.y/2.),
			       (e_outer.z + ef_delta.z/2. + eb_delta.z/2.) };

/* To calculate the alpha decay angle with respect to a vector normal to 
   the plane, we need the normal vector direction cosines.  This is easy, 
   because the detectors are at 90 degrees to each other since they are 
   cosines, we don't care which side of the plane the normal vector sticks out. */

double edgelength=6.4; // Length of the active area edges

double tmx=(w_outer.x-w_later.x)/edgelength; /* top plane normal vector direction 
						cosine wrt x axis */
double tmy=(w_outer.y-w_later.y)/edgelength; /* top plane normal vector direction
						cosine wrt y axis */
double tmz=(w_outer.z-w_later.z)/edgelength; /* top plane normal vector direction 
						cosine wrt z axis */
double emx=(t_outer.x-t_later.x)/edgelength; 
double emy=(t_outer.y-t_later.y)/edgelength; 
double emz=(t_outer.z-t_later.z)/edgelength; 
double wmx=(e_outer.x-e_later.x)/edgelength; 
double wmy=(e_outer.y-e_later.y)/edgelength; 
double wmz=(e_outer.z-e_later.z)/edgelength; 

double umx01=cos(PI*8./6.); /* 1 o'clock upstream plane normal vector 
			       direction cosine wrt x axis */
double umy01=cos(PI*5./6.);
double umz01=0.;	   
double umx03=cos(PI*-6./6.); /* 3 o'clock upstream plane normal vector 
			       direction cosine wrt x axis */
double umy03=cos(PI*-3./6.);
double umz03=0.;		
double umx05=cos(PI* 4./6.);
double umy05=cos(PI* 1./6.);
double umz05=0.;		
double umx07=cos(PI*2./6.);
double umy07=cos(PI*-1./6.);
double umz07=0.;		
double umx09=cos(PI* 0./6.);
double umy09=cos(PI*-3./6.);
double umz09=0.;		
double umx11=cos(PI*-2./6.);
double umy11=cos(PI*-5./6.);
double umz11=0.;		

/* Might as well declare the alpha ray direction cosines here as well... */
double amx,amy,amz;	// direction cosines for alpha path
double deadlayer=0.6;	// normal depth of dead layer in micrometers
double ORdeadlength;	// path length through origin detector dead layer
double TEdeadlength;	// path length through terminus detector dead layer
double ALpathlength;	// distance between origin and terminus

struct rect_point xyz(int chip, int sf, int sb, int ran) {
  /* Takes chip number (0=top, 1=west and 2=east) and strip number from 
     front and back, and returns the (x,y,z) position in centimeters in a 
     rect_point structure.  Bad input returns (0,0,0) indicating a bad
     position.  
     Structure type rect_point must be defined from the scope where this 
     is called.  If ran == 0, the center of the pixel is returned.
     If ran != 0, the positions are randomized over the pixel area. */
  struct rect_point pos = {0.0, 0.0, 0.0};
  double rand0 = 0.0;
  double rand1 = 0.0;

  if (ran) {
    rand0 = (double)rand()/(double)RAND_MAX-0.5;
    rand1 = (double)rand()/(double)RAND_MAX-0.5;
  }

  if ( (chip==0) && (sf>-1) && (sf<nfp) && (sb>-1) && (sb<nfp) ) { // Top
    pos.x = (t_origin.x + ((double)sf + rand0)*tf_delta.x + 
	     ((double)sb + rand1)*tb_delta.x);
    pos.y = (t_origin.y + ((double)sf + rand0)*tf_delta.y + 
	     ((double)sb + rand1)*tb_delta.y);
    pos.z = (t_origin.z + ((double)sf + rand0)*tf_delta.z + 
	     ((double)sb + rand1)*tb_delta.z);
  }
  else if ( (chip==1) && (sf>-1) && (sf<nfp) && (sb>-1) && (sb<nfp) ) { // West
    pos.x = (w_origin.x + ((double)sf + rand0)*wf_delta.x + 
	     ((double)sb + rand1)*wb_delta.x);
    pos.y = (w_origin.y + ((double)sf + rand0)*wf_delta.y + 
	     ((double)sb + rand1)*wb_delta.y);
    pos.z = (w_origin.z + ((double)sf + rand0)*wf_delta.z + 
	     ((double)sb + rand1)*wb_delta.z);
  }
  else if ( (chip==2) && (sf>-1) && (sf<nfp) && (sb>-1) && (sb<nfp) ) { // East
    pos.x = (e_origin.x + ((double)sf + rand0)*ef_delta.x + 
	     ((double)sb + rand1)*eb_delta.x);
    pos.y = (e_origin.y + ((double)sf + rand0)*ef_delta.y + 
	     ((double)sb + rand1)*eb_delta.y);
    pos.z = (e_origin.z + ((double)sf + rand0)*ef_delta.z + 
	     ((double)sb + rand1)*eb_delta.z);
  }
  return pos;
}

#endif
