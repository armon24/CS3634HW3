#include "simpleRayTracer.h"

// HW02Q2: complete this function
bool intersectRayDisk(const ray_t r, const disk_t disk, dfloat *pt_t){

  vector_t s = r.start;
  vector_t d = r.dir;
  vector_t n = disk.normal;
  vector_t c = disk.center;
 
  // A. find if there is an intersection between ray and disk
  //    i.e. exists point p on disk such that p = s+t*d for some t>=0

  // B. If no intersection return false
 
  // C. If there is an intersection then check is t< (*pt_t)
  //    a. if false: return false
  //    b. if true: set *pt_t = t, and return true
 
  double tol = 1e-8;

  //First make sure we intersect the plane containing our disk
  double numerator = vectorDot(n,c) - vectorDot(n,s);

  double denominator = vectorDot(n,d);
 
  if (fabs(denominator) < tol) {
    return false;
  }

  if (numerator/denominator < 0) {
    return false;
  }
 
  //Find the disk/plane intersection
  vector_t test_pt = vectorAdd(s, vectorScale(numerator/denominator, d));

  if (vectorNorm(vectorSub(test_pt, c)) > disk.radius) return false;
 
  double t = numerator/denominator;
 
  if(t < (*pt_t)){
    *pt_t = t;
    return true;
  }

  return false;
}

// HW02Q3: complete this function
bool intersectRayCylinder(const ray_t r, const cylinder_t cylinder, dfloat *pt_t){

  // A. find if there is an intersection between ray and cylinder
  //    i.e. exists point p on disk such that p = s+t*d for some t>=0

  // B. If no intersection return false
 
  // C. If there is an intersection then check is t< (*pt_t)
  //    a. if false: return false
  //    b. if true: set *pt_t = t, and return true

  vector_t s = r.start;
  vector_t d = r.dir;

  vector_t c = cylinder.center;
  vector_t a = cylinder.axis;
  dfloat   R = cylinder.radius;
  dfloat   H = cylinder.height;

  dfloat t = 1e9;
  double ray_dist;

  double tol = 1e-8; //tolerance for floating point numbers
 
  vector_t h_term;               //linear term inside our dot product
  vector_t const_term;           //constant term inside our dot product
  vector_t projected_height;     //vector from the base of the cylinder to the intersection point
  vector_t test_pt;              //guess for the intersection point
 
  //First make sure we intersect an infinite cylinder
  //Use a projection operation to divide a quadratic equation
 
  h_term = vectorSub(d, vectorScale(vectorDot(d,a),a));
            
  const_term = vectorSub(vectorSub(s, c), vectorScale(vectorDot(a, vectorSub(s,c)),a));

  double alpha = vectorDot(h_term,h_term);

  double beta  = 2*vectorDot(h_term,const_term);

  double charlie  = vectorDot(const_term,const_term) - R*R;

  //no intersection at all
  if (beta*beta - 4*alpha*charlie < 0.)
    {
      return false;
    }
 
  //the cylinder is behind us
  if (-1*beta + sqrt(fabs(beta*beta - 4*alpha*charlie)) < tol)
    {
      return false;
    }
 
 
  if (-1*beta > sqrt(beta*beta - 4*alpha*charlie) )
    {
       ray_dist = (-1*sqrt(beta*beta - 4*alpha*charlie) - beta)/(2*alpha);
   
       //test_pt = scaled_sum(1,s,ray_dist,d);
       test_pt = vectorAdd(vectorScale(1, s), vectorScale(ray_dist, d));
   
       //sum of scalar product instead of scaled_sum
       //a + bc (a and c are vectors, b is a scalar
   
       projected_height = vectorSub(test_pt, c);
   
       //check to make sure we are in the bounded cylinder
       if (vectorDot(projected_height,a) > 0 && vectorDot(projected_height,a) < H)
       {
     if(t<(*pt_t))
      {
       *pt_t = t;
       return true;
         }
       }
    }

  ray_dist = (sqrt(beta*beta - 4*alpha*charlie) - beta)/(2*alpha);
 
  //test_pt = scaled_sum(1,s, ray_dist,d);

  test_pt = vectorAdd(vectorScale(1,s), vectorScale(ray_dist, d));
 
  projected_height = vectorSub(test_pt,c);
 
  //check to make sure we are in the bounded cylinder
  if (vectorDot(projected_height,a) > 0 && vectorDot(projected_height,a) < H)
    {
      double t = ray_dist;
     
       if(t< (*pt_t))
        {
         *pt_t = t;
         return true;
        }
    }
 
  return false;

}

// HW02Q4: complete this function
bool intersectRayCone(const ray_t r, const cone_t cone, dfloat *pt_t){

  vector_t s = r.start;
  vector_t d = r.dir;

  vector_t v = cone.vertex;
  vector_t a = cone.axis;
  dfloat   R = cone.radius;
  dfloat   H = cone.height;
 
  dfloat t = 1e9;
 
  // A. find if there is an intersection between ray and cylinder
  //    i.e. exists point p on disk such that p = s+t*d for some t>=0

  // B. If no intersection return false
 
  // C. If there is an intersection then check is t< (*pt_t)
  //    a. if false: return false
  //    b. if true: set *pt_t = t, and return true

  if(t< (*pt_t)){
    *pt_t = t;
    return true;
  }

  return false;
 
}

// Do not edit beyond here--------------------------------------------------------------------------->


/* Check if the ray and triangle intersect */
bool intersectRayTriangle(const ray_t r, const triangle_t tri, dfloat *t){

  // TW: unused fudge factor
  dfloat delta = 0;
 
  bool retval = false;

  vector_t B1 = vectorSub(tri.vertices[2], tri.vertices[0]);
  vector_t B2 = vectorSub(tri.vertices[2], tri.vertices[1]);
  vector_t B3 = r.dir;

  vector_t R = vectorSub(tri.vertices[2], r.start);

  dfloat J = vectorTripleProduct(B2, B3, B1);
 
  dfloat L1 = vectorTripleProduct(B2, B3, R);
  if(L1<delta*J) return false;

  dfloat L2 = vectorTripleProduct(B3, B1, R);
  if(L2<delta*J || L1+L2>J*(1+delta)) return false;

  dfloat t0 = vectorTripleProduct(B1, B2, R)/J;

  /* Verify t1 larger than 0 and less than the original t */
  // TW: FUDGE FACTOR
  if((t0 > p_intersectDelta) && (t0 < *t)){
    *t = t0;
    retval = true;
  }

  return retval;
}

/* Check if the ray and triangle intersect */
bool intersectRayRectangle(const ray_t r, const rectangle_t rect, dfloat *t){

  vector_t C  = rect.center;
  vector_t A1 = rect.axis[0];
  vector_t A2 = rect.axis[1];
  dfloat   L1 = rect.length[0];
  dfloat   L2 = rect.length[1];

  // n = A1 x A2
  // (s + t*d - C).n  = 0
  // t = (C - s).n/(d.n)

  vector_t n = vectorCrossProduct(A1, A2);
  dfloat  t0 = vectorDot(vectorSub(C,r.start), n)/vectorDot(r.dir, n);

  // intersection behind start of ray
  if(t0<0 || t0>*t) return false;

  // X = s + t*d - C
  vector_t X = vectorAdd(vectorSub(r.start,C), vectorScale(t0, r.dir));
 
  dfloat h1 = vectorDot(A1, X)+0.5*L1; // shift
  if(h1<0 || h1>L1) return false;
 
  dfloat h2 = vectorDot(A2, X)+0.5*L2; // shift
  if(h2<0 || h2>L1) return false;

  // success
  *t = t0;
 
  return true;
}



/* Check if the ray and sphere intersect */
bool intersectRaySphere(const ray_t r, const sphere_t s, dfloat *t){
   
  bool retval = false;

  /* A = d.d, the vector_t dot product of the direction */
  dfloat A = vectorDot(r.dir, r.dir);
   
  /* We need a vector_t representing the distance between the start of
   * the ray and the position of the circle.
   * This is the term (p0 - c)
   */
  vector_t dist = vectorSub(r.start, s.pos);
   
  /* 2d.(p0 - c) */ 
  dfloat B = 2.f * vectorDot(r.dir, dist);
   
  /* (p0 - c).(p0 - c) - r^2 */
  dfloat C = vectorDot(dist, dist) - (s.radius * s.radius);

  /* find roots of quadratic */
  dfloat t0, t1;
  if(solveQuadratic(A,0.5*B,C,&t0,&t1)){
    if((t0 > p_intersectDelta) && (t0 < *t)){
      *t = t0;
      retval = true;
    }else
      retval = false;
  }else{
    retval = false;
  }

  return retval;
}


bool intersectPointGridCell(const grid_t   grid,
                const vector_t p,
                const int cellI,
                const int cellJ,
                const int cellK){
 
  if(p.x<=grid.xmin+(cellI  )*grid.dx) return false;
  if(p.x> grid.xmin+(cellI+1)*grid.dx) return false;

  if(p.y<=grid.ymin+(cellJ  )*grid.dy) return false;
  if(p.y> grid.ymin+(cellJ+1)*grid.dy) return false;
 
  if(p.z<=grid.zmin+(cellK  )*grid.dz) return false;
  if(p.z> grid.zmin+(cellK+1)*grid.dz) return false; 
 
  return true;
}

unsigned int intersectRayBox(ray_t *r, const bbox_t bbox){

  vector_t d = r->dir;
  vector_t s = r->start;

  dfloat mint = 20000;
  unsigned int face = 0;
 
  if(d.x>0){ // face 2
    dfloat newt = (bbox.xmax-s.x)/d.x; // d.x > 0
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.x<0){ // face 4
    // s.x + newt*d.x = bbox.xmin
    dfloat newt = (bbox.xmin-s.x)/d.x;
    if(newt>0){
      mint = min(mint, newt);
    }
  }
 
  if(d.y>0){ // face 3
    dfloat newt = (bbox.ymax-s.y)/d.y;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.y<0){ // face 1
    dfloat newt = (bbox.ymin-s.y)/d.y;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.z>0){ // face 5
    dfloat newt = (bbox.zmax-s.z)/d.z;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  if(d.z<0){ // face 0
    dfloat newt = (bbox.zmin-s.z)/d.z;
    if(newt>0){
      mint = min(mint, newt);
    }
  }

  // now figure out which faces the ray passes through
  if(d.x>0){ // face 2
    dfloat newt = (bbox.xmax-s.x)/d.x;
    if(newt>0 && newt<=mint)
      face |= 4;
  }

  if(d.x<0){ // face 4
    dfloat newt = (bbox.xmin-s.x)/d.x;
    if(newt>0 && newt<=mint)
      face |= 16;
  }

  if(d.y>0){ // face 3
    dfloat newt = (bbox.ymax-s.y)/d.y;
    if(newt>0 && newt<=mint)
      face |= 8;
  }

  if(d.y<0){ // face 1
    dfloat newt = (bbox.ymin-s.y)/d.y;
    if(newt>0 && newt<=mint)
      face |= 2;
  }

  if(d.z>0){ // face 5
    dfloat newt = (bbox.zmax-s.z)/d.z;
    if(newt>0 && newt<=mint)
      face |= 32;
  }

  if(d.z<0){ // face 0
    dfloat newt = (bbox.zmin-s.z)/d.z;
    if(newt>0 && newt<=mint)
      face |= 1;
  }
 
  if(face>0){
    // update start of ray
    r->start = vectorAdd(s, vectorScale(mint, d));
  }

  return face;
}

bool intersectRayShape(const ray_t r, const shape_t s, dfloat *t){

  switch(s.type){
  case SPHERE:   return intersectRaySphere  (r, s.sphere,   t);
  case CONE:     return intersectRayCone    (r, s.cone,     t);
  case DISK:     return intersectRayDisk    (r, s.disk,     t);
  case CYLINDER: return intersectRayCylinder(r, s.cylinder, t);
  case RECTANGLE:return intersectRayRectangle(r, s.rectangle, t);
  case TRIANGLE: return intersectRayTriangle(r, s.triangle, t);
  }

  return false;
 
}
