#include <bits/stdc++.h>
#define LL long long
#define endl "\n"
#define EPS 1e-9
#define PI acos(-1.0)
using namespace std;

//=========================0-D data structures ================//
// Point Declaration
//struct Point
//{
//    int x,y;
//};

struct Point_i
{
    int x, y; // whenever possible, work with Point_i
    Point_i()
    {
        x = y = 0;    // default constructor
    }
    Point_i(int _x, int _y) : x(_x), y(_y) {}
}; // user-defined


// Sorting Points
struct Point
{
    double x, y;
    Point()
    {
        x = y = 0.0;
    }
    Point(double _x, double _y) : x(_x), y(_y) {}
    bool operator < (Point other) const   // override less than operator
    {
        if (fabs(x - other.x) > EPS) // useful for sorting
            return x < other.x; // first criteria , by x-coordinate
        return y < other.y;
    }
}; // second criteria, by y-coordinate
//// in int main(), assuming we already have a populated vector<Point> P
//sort(P.begin(), P.end()); // comparison operator is defined above

// Test two Points are equal
//struct Point
//{
//    double x, y;
//    Point()
//    {
//        x = y = 0.0;
//    }
//    Point(double _x, double _y) : x(_x), y(_y) {}
//// use EPS (1e-9) when testing equality of two floating Points
//    bool operator == (Point other) const
//    {
//        return (fabs(x - other.x) < EPS && (fabs(y - other.y) < EPS));
//    }
//};
//// in int main()
//Point P1(0, 0), P2(0, 0), P3(0, 1);
//printf("%d\n", P1 == P2); // true
//printf("%d\n", P1 == P3); // false

// Euclidean distance between Points

double dist(Point p1, Point p2)   // Euclidean distance
{
// hypot(dx, dy) returns sqrt(dx * dx + dy * dy)
    return hypot(p1.x - p2.x, p1.y - p2.y);
} // return double

//rotate a Point by angle8 θ counter clockwise around origin (0, 0) by using a
//rotation matrix:
//
//    |  x' |     |  cos(theta) -sin(theta)  |       |  x  |
//    |     |  =  |                          |   X   |     |
//    | y'  |     |   sin(theta)  cos(theta) |       |  y  |


double DEG_to_RAD(double theta)
{
    return theta*(PI/180.0);
}
double RAD_to_DEG(double theta)
{
    return theta*(180.0/PI);
}


// rotate p by theta degrees CCW w.r.t origin (0, 0)
Point rotate(Point p, double theta)
{
    double rad = DEG_to_RAD(theta); // multiply theta with PI / 180.0
    return Point(p.x * cos(rad) - p.y * sin(rad),
                 p.x * sin(rad) + p.y * cos(rad));
}


//===================1D Objects: Lines =======================

struct line
{
    double a, b, c;
}; // a way to represent a line
// if b==1 --> vertical if b==0 non-vertical

//required line equation

// the answer is stored in the third parameter (pass by reference)
void PointsToLine(Point p1, Point p2, line &l)
{
    if (fabs(p1.x - p2.x) < EPS)   // vertical line is fine
        {
            l.a = 1.0;
            l.b = 0.0;
            l.c = -p1.x; // default values
        }
    else
        {
            l.a = -(double)(p1.y - p2.y) / (p1.x - p2.x);
            l.b = 1.0; // IMPORTANT: we fix the value of b to 1.0
            l.c = -(double)(l.a * p1.x) - p1.y;
        }
}


// parallel and same lines
bool areParallel(line l1, line l2)   // check coefficients a & b
{
    return (fabs(l1.a-l2.a) < EPS) && (fabs(l1.b-l2.b) < EPS);
}
bool areSame(line l1, line l2)   // also check coefficient c
{
    return areParallel(l1 ,l2) && (fabs(l1.c - l2.c) < EPS);
}
// Intersection Point of two lines

// returns true (+ intersection Point) if two lines are intersect
bool areIntersect(line l1, line l2, Point &p)
{
    if (areParallel(l1, l2)) return false; // no intersection
// solve system of 2 linear algebraic equations with 2 unknowns
    p.x = (l2.b * l1.c - l1.b * l2.c) / (l2.a * l1.b - l1.a * l2.b);
// special case: test for vertical line to avoid division by zero
    if (fabs(l1.b) > EPS) p.y = -(l1.a * p.x + l1.c);
    else p.y = -(l2.a * p.x + l2.c);
    return true;
}



// Vectors:
struct vec
{
    double x, y; // name: ‘vec’ is different from STL vector
    vec(double _x, double _y) : x(_x), y(_y) {}
};


vec toVector(Point a, Point b)   // convert 2 Points to vector a->b
{
    return vec(b.x - a.x, b.y - a.y);
}
vec scale(vec v, double s)   // nonnegative s = [<1 .. 1 .. >1]
{
    return vec(v.x * s, v.y * s);
} // shorter.same.longer
Point translate(Point p, vec v)   // translate p according to v
{
    return Point(p.x + v.x , p.y + v.y);
}


/*
Given a Point p and a line l (described by two Points a and b), we can compute the
minimum distance from p to l by first computing the location of Point c in l that is
closest to Point p (see Figure 7.2—left) and then obtain the Euclidean distance between
p and c.
*/

double dot(vec a, vec b)
{
    return (a.x * b.x + a.y * b.y);
}
double norm_sq(vec v)
{
    return v.x * v.x + v.y * v.y;
}
// returns the distance from p to the line defined by
// two Points a and b (a and b must be different)
// the closest Point is stored in the 4th parameter (byref)
double distToLine(Point p, Point a, Point b, Point &c)
{
// formula: c = a + u * ab
    vec ap = toVector(a, p), ab = toVector(a, b);
    double u = dot(ap, ab) / norm_sq(ab);
    c = translate(a, scale(ab, u)); // translate a to c
    return dist(p, c);
} // Euclidean distance between p and c

/*If we are given a line segment instead (defined by two end Points a and b), then the
minimum distance from Point p to line segment ab must also consider two special
cases, the end Points a and b of that line segment*/

// returns the distance from p to the line segment ab defined by
// two Points a and b (still OK if a == b)
// the closest Point is stored in the 4th parameter (byref)
double distToLineSegment(Point p, Point a, Point b, Point &c)
{
    vec ap = toVector(a, p), ab = toVector(a, b);
    double u = dot(ap, ab) / norm_sq(ab);
    if (u < 0.0)
        {
            c = Point(a.x, a.y); // closer to a
            return dist(p, a);
        } // Euclidean distance between p and a
    if (u > 1.0)
        {
            c = Point(b.x, b.y); // closer to b
            return dist(p, b);
        } // Euclidean distance between p and b
    return distToLine(p, a, b, c);
} // run distToLine as above


/* the angle aob given three Points: a, o, and b, using dot product13.
Since oa · ob = |oa| × |ob| × cos(θ), we have theta = arccos(oa · ob/(|oa| × |ob|)).
*/

double angle(Point a, Point o, Point b)   // returns angle aob in rad
{
    vec oa = toVector(o, a);
    vec ob = toVector(o, b);
    return acos(dot(oa, ob) / sqrt(norm_sq(oa) * norm_sq(ob)));
}

/* Given a line defined by two Points p and q, we can determine whether a Point r is
on the left/right side of the line, or whether the three Points p, q, and r are collinear.
This can be determined with cross product. Let pq and pr be the two vectors obtained
from these three Points. The cross product pq × pr result in another vector that is
perpendicular to both pq and pr. The magnitude of this vector is equal to the area of
the parallelogram that the vectors span14. If the magnitude is positive/zero/negative,
then we know that p → q → r is a left turn/collinear/right turn, respectively.
 The left turn test is more famously known as the CCW (Counter
Clockwise) Test.*/


double cross(vec a, vec b)
{
    return a.x * b.y - a.y * b.x;
}
// note: to accept collinear Points, we have to change the '> 0'
// returns true if Point r is on the left side of line pq
bool ccw(Point p, Point q, Point r)
{
    return cross(toVector(p, q), toVector(p, r)) > 0;
}
// returns true if Point r is on the same line as the line pq
bool collinear(Point p, Point q, Point r)
{
    return fabs(cross(toVector(p, q), toVector(p, r))) < EPS;
}

//============== Line using slope and intercept form===========
struct line2
{
    double m, c;
}; // another way to represent a line
int PointsToLine2(Point p1, Point p2, line2 &l)
{
    if (p1.x == p2.x)   // special case: vertical line
        {
            l.m = INF; // l contains m = INF and c = x_value
            l.c = p1.x; // to denote vertical line x = x_value
            return 0; // we need this return variable to differentiate result
        }
    else
        {
            l.m = (double)(p1.y - p2.y) / (p1.x - p2.x);
            l.c = p1.y - l.m * p1.x;
            return 1; // l contains m and c of the line equation y = mx + c
        }
}

//to compute line equation given a Point and gradient
// convert Point and gradient/slope to line
void PointSlopeToLine(Point p, double m, line &l)
{
    l.a = -m; // always -m
    l.b = 1; // always 1
    l.c = -((l.a * p.x) + (l.b * p.y));
} // compute this


void closestPoint(line l, Point p, Point &ans)
{
    line perpendicular; // perpendicular to l and pass through p
    if (fabs(l.b) < EPS)   // special case 1: vertical line
        {
            ans.x = -(l.c);
            ans.y = p.y;
            return;
        }
    if (fabs(l.a) < EPS)   // special case 2: horizontal line
        {
            ans.x = p.x;
            ans.y = -(l.c);
            return;
        }
    PointSlopeToLine(p, 1 / l.a, perpendicular); // normal line
// intersect line l with this perpendicular line
// the intersection Point is the closest Point
    areIntersect(l, perpendicular, ans);
}


// returns the reflection of Point on a line
void reflectionPoint(line l, Point p, Point &ans)
{
    Point b;
    closestPoint(l, p, b); // similar to distToLine
    vec v = toVector(p, b); // create a vector
    ans = translate(translate(p, v), v);
} // translate p twice



// ================ 2 D Object circles ===============

// To check whether a Point p is inside a circle c of radius r or not

int insideCircle(Point_i p, Point_i c, int r)   // all integer version
{
    int dx = p.x - c.x, dy = p.y - c.y;
    int Euc = dx * dx + dy * dy, rSq = r * r; // all integer
    return Euc < rSq ? 0 : Euc == rSq ? 1 : 2;
} //inside/border/outside

// Given a  central angle  theta , and circumference ci = 2*pi*r
// Length Arc of a circle = (2*pi*r * theta)/360.0;
// Length of chord =  sqrt(2*r*r*(1- cos(theta))) OR  2*r*sin(theta/2)))
// Area of circle = pi*r*r
// Sector area with central angle theta = (theta*area)/360.0;
/*Segment of a circle is defined as a region of the circle enclosed by a chord and an arc
lying between the chord’s endPoints . The area of a segment can be found by subtracting the area of the corresponding sector from the area of an
isosceles triangle with sides: r, r, and chord-length.*/


/*Given 2 Points on the circle (p1 and p2) and radius r of the corresponding circle, we
can determine the location of the centers (c1 and c2) of the two possible circles */
bool circle2PtsRad(Point p1, Point p2, double r, Point &c)
{
    double d2 = (p1.x - p2.x) * (p1.x - p2.x) +
                (p1.y - p2.y) * (p1.y - p2.y);
    double det = r * r / d2 - 0.25;
    if (det < 0.0) return false;
    double h = sqrt(det);
    c.x = (p1.x + p2.x) * 0.5 + (p1.y - p2.y) * h;
    c.y = (p1.y + p2.y) * 0.5 + (p2.x - p1.x) * h;
    return true;
} // to get the other center, reverse p1 and p2


//=================2D Objects: Triangles =======================

double perimeter(double ab,double bc,double ca)
{
    return ab+bc+ca;
}
double perimeter(Point a, Point b, Point c)
{
    return dist(a, b) + dist(b, c) + dist(c, a);
}

double area(double ab,double bc,double ca)
{
    double s = perimeter(ab,bc,ca)/2;
    return sqrt(s*(s-ab)*(s-bc)*(s-ca));
}
double area(Point a, Point b, Point c)
{
    return area(dist(a, b), dist(b, c), dist(c, a));
}

/*A triangle with area A and semi-perimeter s has an inscribed circle (incircle) with
radius r = A/s */
double rInCircle(double ab, double bc, double ca)
{
    return area(ab, bc, ca) / (0.5 * perimeter(ab, bc, ca));
}
double rInCircle(Point a, Point b, Point c)
{
    return rInCircle(dist(a, b), dist(b, c), dist(c, a));
}

/* Incircle The center of incircle is the meeting Point between the triangle's angle bisectors . We can get the center if we have two angle bisectors and find their
intersection Point.  */

// assumption: the required Points/lines functions have been written
// returns 1 if there is an inCircle center, returns 0 otherwise
// if this function returns 1, ctr will be the inCircle center
// and r is the same as rInCircle
int inCircle(Point p1, Point p2, Point p3, Point &ctr, double &r)
{
    r = rInCircle(p1, p2, p3);
    if (fabs(r) < EPS) return 0; // no inCircle center
    line l1, l2; // compute these two angle bisectors
    double ratio = dist(p1, p2) / dist(p1, p3);
    Point p = translate(p2, scale(toVector(p2, p3), ratio / (1 + ratio)));
    PointsToLine(p1, p, l1);
    ratio = dist(p2, p1) / dist(p2, p3);
    p = translate(p1, scale(toVector(p1, p3), ratio / (1 + ratio)));
    PointsToLine(p2, p, l2);
    areIntersect(l1, l2, ctr); // get their intersection Point
    return 1;
}



/*A triangle with 3 sides: a, b, c and area A has an circumscribed circle
(circumcircle) with radius R = a × b × c/(4 × A)*/

double rCircumCircle(double ab, double bc, double ca)
{
    return ab * bc * ca / (4.0 * area(ab, bc, ca));
}
double rCircumCircle(Point a, Point b, Point c)
{
    return rCircumCircle(dist(a, b), dist(b, c), dist(c, a));
}

bool canFormTriangle(double a, double b, double c)
{
    return (a + b > c) && (a + c > b) && (b + c > a);
}

//===================Algorithm on Polygon with Libraries==============================//

/*perimeter of a polygon (either convex or concave) with n vertices given in some order
(either clockwise or counter-clockwise)*/

// returns the perimeter, which is the sum of Euclidian distances
// of consecutive line segments (polygon edges)
double perimeter(const vector<Point> &P)
{
    double result = 0.0;
    for (int i = 0; i < (int)P.size()-1; i++) // remember that P[0] = P[n-1]
        result += dist(P[i], P[i+1]);
    return result;
}


// Area of Polygon
// returns the area, which is half the determinant
double area(const vector<Point> &P)
{
    double result = 0.0, x1, y1, x2, y2;
    for (int i = 0; i < (int)P.size()-1; i++)
        {
            x1 = P[i].x;
            x2 = P[i+1].x;
            y1 = P[i].y;
            y2 = P[i+1].y;
            result += (x1 * y2 - x2 * y1);
        }
    return fabs(result) / 2.0;
}

// Checking whether a polygon is convex or concave
bool isConvex(const vector<Point> &P)   // returns true if all three
{
    int sz = (int)P.size(); // consecutive vertices of P form the same turns
    if (sz <= 3) return false; // a Point/sz=2 or a line/sz=3 is not convex
    bool isLeft = ccw(P[0], P[1], P[2]); // remember one result
    for (int i = 1; i < sz-1; i++) // then compare with the others
        if (ccw(P[i], P[i+1], P[(i+2) == sz ? 1 : i+2]) != isLeft)
            return false; // different sign -> this polygon is concave
    return true; // this polygon is convex
}


// Checking whether a Point is inside a 2D polygon

// wn_PnPoly(): winding number test for a Point in a polygon
//      Input:   P = a Point,
//               V[] = vertex Points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int wn_PnPoly( Point P, const vector<Point> &V, int n )
{
    int    wn = 0;    // the  winding number counter
    // loop through all edges of the polygon
    for (int i=0; i<n; i++)     // edge from V[i] to  V[i+1]
        {
            if (V[i].y <= P.y)            // start y <= P.y
                {
                    if (V[i+1].y  > P.y)      // an upward crossing
                        if (ccw( V[i], V[i+1], P)==true)  // P left of  edge
                            ++wn;            // have  a valid up intersect
                }
            else                          // start y > P.y (no test needed)
                {
                    if (V[i+1].y  <= P.y)     // a downward crossing
                        if (ccw( V[i], V[i+1], P) ==false)  // P right of  edge
                            --wn;            // have  a valid down intersect
                }
        }
    return wn;
}
//===================Cutting Polygon with a straight line=====================
// line segment p-q intersect with line A-B.
Point lineIntersectSeg(Point p, Point q, Point A, Point B)
{
    double a = B.y - A.y;
    double b = A.x - B.x;
    double c = B.x * A.y - A.x * B.y;
    double u = fabs(a * p.x + b * p.y + c);
    double v = fabs(a * q.x + b * q.y + c);
    return Point((p.x * v + q.x * u) / (u+v), (p.y * v + q.y * u) / (u+v));
}

// cuts polygon Q along the line formed by Point a -> Point b
// (note: the last Point must be the same as the first Point)
vector<Point> cutPolygon(Point a, Point b, const vector<Point> &Q)
{
    vector<Point> P;
    for (int i = 0; i < (int)Q.size(); i++)
        {
            double left1 = cross(toVector(a,b),toVector(a,Q[i])), left2 = 0;
            if (i != (int)Q.size()-1)
                left2 = cross(toVector(a,b),toVector(a,Q[i+1]));
            if (left1 > -EPS) P.push_back(Q[i]); // Q[i] is on the left of ab
            if (left1 * left2 < -EPS) // edge (Q[i], Q[i+1]) crosses line ab
                P.push_back(lineIntersectSeg(Q[i], Q[i+1], a, b));
        }
    if (!P.empty() && (fabs(P.back().x - P.front().x) > EPS ||
                       fabs(P.back().y - P.front().y) > EPS))
        P.push_back(P.front()); // make P's first Point = P's last Point
    return P;
}
//===================Convex Hull======================
Point pivot(0, 0);
bool angle_cmp(Point a, Point b)   // angle-sorting function
{
    if (collinear(pivot, a, b))
        return dist(pivot, a) < dist(pivot, b); // which one is closer?
    double d1x = a.x - pivot.x, d1y = a.y - pivot.y;
    double d2x = b.x - pivot.x, d2y = b.y - pivot.y;
    return (atan2(d1y, d1x) - atan2(d2y, d2x)) < 0;
}

vector<Point> CH(vector<Point> &P)   // P can change after sort
{
    int i, j, n = (int)P.size();
    if (n <= 3) return P; // special case, the CH is P itself
    // first, find P0 = Point with lowest Y and if tie: rightmost X
    int P0 = 0;
    for (i = 1; i < n; i++) // O(n)
        if (P[i].y  < P[P0].y ||
                (P[i].y == P[P0].y && P[i].x > P[P0].x))
            P0 = i;
    // swap selected vertex with P[0]
    Point temp = P[0];
    P[0] = P[P0];
    P[P0] = temp;
    // second, sort Points by angle w.r.t. P0, skipping P[0]
    pivot = P[0]; // use this global variable as reference
    sort(++P.begin(), P.end(), angle_cmp); // O(n log n) here
    // third, the ccw tests
    Point prev(0, 0), now(0, 0);
    vector<Point> S;
    S.push_back(P[n-1]);
    S.push_back(P[0]);
    S.push_back(P[1]); // initial S
    i = 2; // then, we check the rest
    while (i < n)   // note: n must be >= 3 for this method to work
        {
            j = (int)S.size() - 1;
            if (ccw(S[j-1], S[j], P[i])) S.push_back(P[i++]); // left turn, accept
            else S.pop_back(); // or pop the top of S until we have a left turn
        }
    return S;
} // return the result


int main()
{
    // MAIN FOR 0-Dimensional
    Point p1(2.3,2),p2(6,5);
    cout<<" Euclid distance between Points : "<<dist(p1,p2)<<endl;
    Point p3(10,3);
    Point p3_rotated;
    p3_rotated= rotate(p3,90);
    cout<<" Rotated Points : "<<p3_rotated.x<<" , "<<p3_rotated.y<<endl;
    p3_rotated= rotate(p3,77);
    cout<<" Rotated Points : "<<p3_rotated.x<<" , "<<p3_rotated.y<<endl;
    vector<Point> P(3);
    P[0].x = 2;
    P[0].y = 1;
    P[1].x = 4;
    P[1].y = 2;
    P[2].x = 1;
    P[2].y = 4;
    sort(P.begin(),P.end());
    cout<<" Sorted  Points :"<<endl;
    for(vector<Point>::iterator it=P.begin(); it!=P.end(); it++)
        {
            cout<<" "<<(*it).x<<","<<(*it).y<<endl;
        }
    cout<<endl;
    // ============== 1D Main =============
    // Exercise 7.2.2.2: Compute line equation that pass through two Points (2, 2) and (4, 3)!
    Point p4(2,2),p5(4,3);
    line l1;
    PointsToLine(p4,p5,l1);
    cout<<" Coefficients of line l1 : "<<endl;
    cout<<" a : "<<l1.a<< " b : "<<l1.b<<" c : "<<l1.c<<endl;
    /*Exercise 7.2.2.6: Translate a Point c (3, 2) according to a vector ab defined by two Points:
    a (2, 2) and b (4, 3). What is the new coordinate of the Point?*/
    Point c2(3,2),a1(2,2),b1(4,3);
    vec ab = toVector(a1,b1);
    Point ans = rotate(translate(c2,scale(ab,0.5)),90);
    cout<<" New Coordinate : "<<ans.x<<" , "<<ans.y<<endl;
// ================== 2 D Main =========================
    Point_i p6(2,3),p7(-1,5),p8(10,4),c1(4,5);
    cout<<" insideCircle(p6,c1,5) "<<insideCircle(p6,c1,5)<<endl;
    cout<<" insideCircle(p7,c1,5) "<<insideCircle(p7,c1,5)<<endl;
    cout<<" insideCircle(p8,c1,5) "<<insideCircle(p8,c1,5)<<endl;
    //=== 2D circles
    // circle equation, inside, border, outside
    Point_i pt(2, 2);
    int r = 7;
    Point_i inside(8, 2);
    printf("%d\n", insideCircle(inside, pt, r)); // 0-inside
    Point_i border(9, 2);
    printf("%d\n", insideCircle(border, pt, r)); // 1-at border
    Point_i outside(10, 2);
    printf("%d\n", insideCircle(outside, pt, r)); // 2-outside
    double d = 2 * r;
    printf("Diameter = %.2lf\n", d);
    double ci = PI * d;
    printf("Circumference / Perimeter = %.2lf\n", ci);
    double A = PI * r * r;
    printf("Area of circle = %.2lf\n", A);
    printf("Length of arc (central angle = 30 degrees) = %.2lf\n", 30.0 / 360.0 * ci);
    printf("Length of chord (central angle = 30 degrees) = %.2lf\n", sqrt((2 * r * r) * (1 - cos(DEG_to_RAD(30.0)))));
    printf("Area of sector (central angle = 30 degrees) = %.2lf\n", 30.0 / 360.0 * A);
    // printf("Area of segment (30 degrees) = %.2lf\n", 0.0); // subtract the area of the corresponding sector with the area of an Isosceles Triangle with sides: r, r, and Chord-length
    Point p15(0.0,  1.0);
    Point p16(0.0, -1.0);
    Point ans1(0.0, 0.0);
    if(circle2PtsRad(p15, p16, 2.0, ans1))
        printf("One of the center is (%.2lf, %.2lf)\n", ans1.x, ans1.y);
    else cout<<" Not Possible ";
    if(circle2PtsRad(p16, p15, 2.0, ans1))
        printf("The other center  is (%.2lf, %.2lf)\n", ans1.x, ans1.y);
    else cout<<" Not Possible ";
//================Triangles======================
    double base = 4.0, h = 3.0;
    double Ar = 0.5 * base * h;
    printf("Area = %.2lf\n", Ar);
    Point a(0.0, 0.0); // a right triangle
    Point b(4.0, 0.0);
    Point c(4.0, 3.0);
    double p = perimeter(a, b, c);
    double s = 0.5 * p;
    Ar = area(a, b, c);
    printf("Area = %.2lf\n", Ar); // must be the same as above
    printf("R1 (radius of inner circle) = %.2lf\n", rInCircle(a, b, c));
    printf("R2 (radius of outer circle) = %.2lf\n", rCircumCircle(a, b, c));
    // Law of Cosines
    double ab1 = dist(a, b);
    double bc1 = dist(b, c);
    double ca1 = dist(c, a);
    double alpha = RAD_to_DEG(acos((ca1 * ca1 + ab1 * ab1 - bc1 * bc1) / (2.0 * ca1 * ab1)));
    printf("alpha = %.2lf\n", alpha);
    double beta  = RAD_to_DEG(acos((ab1 * ab1 + bc1 * bc1 - ca1 * ca1) / (2.0 * ab1 * bc1)));
    printf("beta  = %.2lf\n", beta);
    double gamma = RAD_to_DEG(acos((bc1 * bc1 + ca1 * ca1 - ab1 * ab1) / (2.0 * bc1 * ca1)));
    printf("gamma = %.2lf\n", gamma);
    // Law of Sines
    printf("%.2lf == %.2lf == %.2lf\n", bc1 / sin(DEG_to_RAD(alpha)), ca1 / sin(DEG_to_RAD(beta)), ab1 / sin(DEG_to_RAD(gamma)));
    // Phytagorean Theorem
    printf("%.2lf^2 == %.2lf^2 + %.2lf^2\n", ca1, ab1, bc1);
    // Triangle Inequality
    printf("(%d, %d, %d) => can form triangle? %d\n", 3, 4, 5, canFormTriangle(3, 4, 5)); // yes
    printf("(%d, %d, %d) => can form triangle? %d\n", 3, 4, 7, canFormTriangle(3, 4, 7)); // no, actually straight line
    printf("(%d, %d, %d) => can form triangle? %d\n", 3, 4, 8, canFormTriangle(3, 4, 8)); // no
// =========================== Algortihms for polygon libraries==================
// 6 Points, entered in counter clockwise order, 0-based indexing
    vector<Point> Pa;
    Pa.push_back(Point(1, 1)); // P0
    Pa.push_back(Point(3, 3)); // P1
    Pa.push_back(Point(9, 1)); // P2
    Pa.push_back(Point(12, 4)); // P3
    Pa.push_back(Point(9, 7)); // P4
    Pa.push_back(Point(1, 7)); // P5
    Pa.push_back(Pa[0]); // important: loop back
    printf("Perimeter of polygon = %.2lf\n", perimeter(Pa)); // 31.64
    printf("Area of polygon = %.2lf\n", area(Pa)); // 49.00
    printf("Is convex = %d\n", isConvex(Pa)); // false (P1 is the culprit)
    //// the positions of P6 and P7 w.r.t the polygon
    //7 P5--------------P4
    //6 |                  \
    //5 |                    \
    //4 |   P7                P3
    //3 |   P1___            /
    //2 | / P6    \ ___    /
    //1 P0              P2
    //0 1 2 3 4 5 6 7 8 9 101112
    Point P6(3, 2); // outside this (concave) polygon
    printf("Point P6 is inside this polygon = %d\n", wn_PnPoly(P6, Pa,Pa.size())); // false =  0 only when P is outside
    Point P7(3, 4); // inside this (concave) polygon
    printf("Point P7 is inside this polygon = %d\n", wn_PnPoly(P7, Pa,Pa.size())); // true = 1 if Point is inside Polygon
    // cutting the original polygon based on line P[2] -> P[4] (get the left side)
    //7 P5--------------P4
    //6 |               |  \
    //5 |               |    \
    //4 |               |     P3
    //3 |   P1___       |    /
    //2 | /       \ ___ |  /
    //1 P0              P2
    //0 1 2 3 4 5 6 7 8 9 101112
    // new polygon (notice the index are different now):
    //7 P4--------------P3
    //6 |               |
    //5 |               |
    //4 |               |
    //3 |   P1___       |
    //2 | /       \ ___ |
    //1 P0              P2
    //0 1 2 3 4 5 6 7 8 9
    Pa = cutPolygon(Pa[2], Pa[4], Pa);
    printf("Perimeter of polygon = %.2lf\n", perimeter(Pa)); // smaller now 29.15
    printf("Area of polygon = %.2lf\n", area(Pa)); // 40.00
    // running convex hull of the resulting polygon (index changes again)
    //7 P3--------------P2
    //6 |               |
    //5 |               |
    //4 |   P7          |
    //3 |               |
    //2 |               |
    //1 P0--------------P1
    //0 1 2 3 4 5 6 7 8 9
    Pa = CH(Pa); // now this is a rectangle
    printf("Perimeter of polygon = %.2lf\n", perimeter(Pa)); // precisely 28.00
    printf("Area of polygon = %.2lf\n", area(Pa)); // precisely 48.00
    printf("Is convex = %d\n", isConvex(Pa)); // true
    printf("Point P7 is inside this polygon = %d\n", wn_PnPoly(P7, Pa,Pa.size())); // true
    return 0;
}
