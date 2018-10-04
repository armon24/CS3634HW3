#include <math.h>
#include "intersections.h"
#include <stdio.h>

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_sphere_intersection(ray_t observer, sphere_t obj, vector_t *intersection) {
  
  //linear term inside the norm
  vector_t h_term = observer.dir;
  //constant term inside the norm
  vector_t const_term = difference(observer.start,obj.center);
  
  //quadratic term in the expansion
  double a = dot_product(h_term,h_term);
  //linear term in the expansion
  double b = 2*dot_product(h_term,const_term);
  //constant term in the expansion
  double c = dot_product(const_term,const_term) - obj.radius*obj.radius;

  double discriminant = b*b - 4*a*c;

  //make sure the equation has real roots
  if (discriminant < 0.) return 0;
  //find the smallest positive h
  double h;
  if ((-1*b - sqrt(discriminant))/(2*a) > 0) h = (-1*b - sqrt(discriminant))/(2*a); //sphere in front of us
  else if ((-1*b + sqrt(discriminant))/(2*a) > 0) h = (-1*b + sqrt(discriminant))/(2*a); //we are inside the sphere
  else return 0; //no positive roots, so the sphere is behind us

  vector_t solution = scaled_sum(1.,observer.start,h,observer.dir);

  copy_vector(solution,intersection);
  return 1;  
}

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_disk_intersection(ray_t observer, disk_t obj, vector_t *intersection) {
  
  //dot_product takes two vector's as parameters. The following dot product tests if the ray runs parallel to the disk
  if(dot_product(observer.dir,obj.normal) == 0)
  {
	return 0;
  }
  
  //the following will test that the array is facing forward and sees the disk. If it does not, that means
  //that the disk is behind the start of the array
 
  double num  = dot_product(obj.normal,difference(obj.center,observer.start)); 
  
  double denom = dot_product(observer.dir, obj.normal);

  double t = (num/denom);

  if(t < 0.0)
  {
	return 0;
  }

  //the following tests the final case--that the array is going through the center of the disk
  if( magnitude( difference( sum( observer.start, scalar_product(t, observer.dir) ), obj.center ) ) < obj.radius )
  {
        //fills the pointer vector with the correct attributes and returns 1, indicating an intersection
	copy_vector( sum( observer.start, scalar_product(t, observer.dir) ) , intersection);
	return 1;
  }
  else
  {
	return 0;
  }
}

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_cylinder_intersection(ray_t observer, cylinder_t obj, vector_t *intersection) {
  //Question 5: Modify this function to compute an intersection
  
  //solve for T 

  //firsthalf of the Projection Operation
  //double firstHalf = difference(sum(observer.start, scalar_product(t, observer.dir)), obj.center); // p = sum(observer.start, scalar_product(t, observer.dir))

 // double secondHalf = magnitude( cross_product( dot_product( obj.axis, difference(sum(observer.start, scalar_product(t, observer.dir)),obj.center) , obj.axis )   ); 

  //cross_product( dot_product( obj.axis, difference(sum(observer.start, scalar_product(t, observer.dir)),obj.center) , obj.axis  )
  
 // if( firstHalf - secondHalf < obj.radius )
 // {
//	return 1;
 // }
 // else
 // {
//	return 0;
 //
  


 //initializes a new vector for the first part of solving for t
 
 vector_t firstParam = difference(observer.start, obj.center) ;
 
//vector_t

 double k = dot_product(obj.axis, difference(observer.start, obj.center));

 vector_t newVec;
 newVec.x = k * obj.axis.x;
 newVec.y = k * obj.axis.y;
 newVec.z = k * obj.axis.z; 

 vector_t secondParam = cross_product(newVec, obj.axis);

 //double  secondParam = magnitude( cross_product(obj.axis, dot_product(obj.axis, difference(observer.start, obj.center))) ); 

 vector_t x = difference(firstParam,secondParam);

 double temp = dot_product(obj.axis, observer.dir);

 vector_t newVec2;
 newVec2.x = temp * obj.axis.x;
 newVec2.y = temp * obj.axis.y;
 newVec2.z = temp * obj.axis.z;

 vector_t firstPart = cross_product(newVec2, obj.axis);

 vector_t y = difference(firstPart, observer.dir);

 //quadratic formula pieces
 double  a = pow( magnitude(y), 2 );

 double  b = (-2) * ( dot_product(x, y) );

 double c = magnitude(x) - pow(obj.radius, 2);
 
 //finally solving for T
 double firstHalfT = (-1)*b;

 double secondHalfT = sqrt( pow(b, 2) - (4*a*c) );

 double t1 = (firstHalfT + secondHalfT) / (2*a);

 double t2 = (firstHalfT - secondHalfT) / (2*a);

 double discriminant = pow(b, 2) - (4*a*c);

 if( discriminant == 0 )
 {
	return 1;
 }
 else
 {
        return 0;
 }

// Steps
//(1) Solve for T
//(2) Solve for A using projection operations
//(3) Do dot product of A and observer.dir
//(4) test if T < 0 and check if ray is facing the cylinder
//(5) check inside? i don't think so
}

int ray_cone_intersection(ray_t observer, cone_t obj, vector_t *intersection) {
  //Question 7: Modify this function to compute an intersection
  return 0;
}
