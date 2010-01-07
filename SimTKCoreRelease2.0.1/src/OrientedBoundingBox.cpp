/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "simbody/internal/OrientedBoundingBox.h"
#include "SimTKmath.h"

namespace SimTK {

OrientedBoundingBox::OrientedBoundingBox() {
}

OrientedBoundingBox::OrientedBoundingBox(const Transform& transform, const Vec3& size) : transform(transform), size(size) {
}

const Transform& OrientedBoundingBox::getTransform() const {
    return transform;
}

const Vec3& OrientedBoundingBox::getSize() const {
    return size;
}

OrientedBoundingBox::OrientedBoundingBox(const Vector_<Vec3>& points) {
    SimTK_APIARGCHECK(points.size() > 0, "OrientedBoundingBox", "OrientedBoundingBox", "No points passed to constructor");
    
    // Construct the covariance matrix of the points.
    
    Vec3 center = mean(points);
    Vector_<Vec3> p = points-center;
    Mat33 c(0);
    for (int i = 0; i < p.size(); i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                c(j, k) += p[i][j]*p[i][k];
    c *= 1.0/p.size();
    
    // Find the eigenvectors, which will be our initial guess for the axes of the box.
    
    Vector_<std::complex<Real> > eigenvalues;
    Matrix_<std::complex<Real> > eigenvectors;
    Eigen(Matrix(c)).getAllEigenValuesAndVectors(eigenvalues, eigenvectors);
    Vec3 axes[3];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
            axes[i][j] = eigenvectors(j, i).real();

    // Now try optimizing the rotation to give a better fit.
    
    Rotation rot(UnitVec3(axes[0]), XAxis, axes[1], YAxis);
    Real volume = calculateVolume(points, rot);
    for (Real step = 0.1; step > 0.01; step *= 0.5) {
        bool improved = true;
        while (improved) {
            Rotation trialRotation[6];
            trialRotation[0].setRotationFromAngleAboutX(step);
            trialRotation[1].setRotationFromAngleAboutX(-step);
            trialRotation[2].setRotationFromAngleAboutY(step);
            trialRotation[3].setRotationFromAngleAboutY(-step);
            trialRotation[4].setRotationFromAngleAboutZ(step);
            trialRotation[5].setRotationFromAngleAboutZ(-step);
            improved = false;
            for (int i = 0; i < 6; i++) {
                trialRotation[i] = trialRotation[i]*rot;
                Real trialVolume = calculateVolume(points, trialRotation[i]);
                if (trialVolume < volume) {
                    rot = trialRotation[i];
                    volume = trialVolume;
                    improved = true;
                }
            }
        }
    }
    
    // Find the extent along each axis.
  
    axes[0] = Vec3(rot.col(0));
    axes[1] = Vec3(rot.col(1));
    axes[2] = Vec3(rot.col(2));
    Vec3 minExtent = Vec3(MostPositiveReal);
    Vec3 maxExtent = Vec3(MostNegativeReal);
    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < 3; j++) {
            minExtent[j] = std::min(minExtent[j], ~axes[j]*points[i]);
            maxExtent[j] = std::max(maxExtent[j], ~axes[j]*points[i]);
        }
    }
    
    // Create the bounding box.
    
    size = maxExtent-minExtent;
    Vec3 tol = 1e-5*size;
    for (int i = 0; i < 3; i++)
        tol[i] = std::max(tol[i], 1e-10);
    size += 2*tol;
    transform = Transform(rot, rot*(minExtent-tol));
}

Real OrientedBoundingBox::calculateVolume(const Vector_<Vec3>& points, const Rotation& rotation) {
    Vec3 minExtent = Vec3(MostPositiveReal);
    Vec3 maxExtent = Vec3(MostNegativeReal);
    for (int i = 0; i < points.size(); i++) {
        Vec3 p = ~rotation*points[i];
        for (int j = 0; j < 3; j++) {
            minExtent[j] = std::min(minExtent[j], p[j]);
            maxExtent[j] = std::max(maxExtent[j], p[j]);
        }
    }
    Vec3 size = maxExtent-minExtent+Vec3(2e-10);
    return size[0]*size[1]*size[2];
}

bool OrientedBoundingBox::containsPoint(const Vec3& point) const {
    Vec3 p = ~transform*point;
    return (p[0] >= 0 && p[0] <= size[0] &&
            p[1] >= 0 && p[1] <= size[1] &&
            p[2] >= 0 && p[2] <= size[2]);
}

bool OrientedBoundingBox::intersectsBox(const OrientedBoundingBox& box) const {
    // Precalculate various quantities.
    
    const Transform t = ~getTransform()*box.getTransform(); // From the other box's frame to this one's
    const Mat33& r = t.R().asMat33();
    const Mat33 rabs = r.abs();
    const Vec3 a = 0.5*getSize();
    const Vec3 b = 0.5*box.getSize();
    const Vec3 center1 = a;
    const Vec3 center2 = t*b;
    const Vec3 d = center2-center1;
    
    // Now perform a series of 15 tests where we project each box onto an axis and see if
    // they overlap.  This is described in Gottschalk, S., Lin, MC, Manocha, D, "OBBTree:
    // a hierarchical structure for rapid interference detection." Proceedings of the 23rd
    // Annual Conference on Computer Graphics and Interactive Techniques, pp. 171-180, 1996.
    // We also perform an additional check which allows an early acceptance if the center of
    // one box is inside the other one.
    
    // First check the three axes of this box.
    
    bool accept = true;
    for (int i = 0; i < 3; i++) {
        Real ra = a[i];
        Real rb = rabs.row(i)*b;
        Real distance = std::abs(d[i]);
        if (distance > ra+rb)
            return false;
        if (distance > ra)
            accept = false;
    }
    if (accept)
        return true;
    
    // Now check the three axes of the other box.
    
    accept = true;
    for (int i = 0; i < 3; i++) {
        Real ra = ~a*rabs.col(i);
        Real rb = b[i];
        Real distance = std::abs(d[0]*r(0, i)+d[1]*r(1, i)+d[2]*r(2, i));
        if (distance > ra+rb)
            return false;
        if (distance > rb)
            accept = false;
    }
    if (accept)
        return true;
    
    // Now check the nine axes formed from cross products of one axis from each box.
    
    {
        Real ra = a[1]*rabs(2, 0)+a[2]*rabs(1, 0);
        Real rb = b[1]*rabs(0, 2)+b[2]*rabs(0, 1);
        if (std::abs(d[2]*r(1, 0) - d[1]*r(2, 0)) > ra+rb)
            return false;
    }
    {
        Real ra = a[1]*rabs(2, 1)+a[2]*rabs(1, 1);
        Real rb = b[0]*rabs(0, 2)+b[2]*rabs(0, 0);
        if (std::abs(d[2]*r(1, 1) - d[1]*r(2, 1)) > ra+rb)
            return false;
    }
    {
        Real ra = a[1]*rabs(2, 2)+a[2]*rabs(1, 2);
        Real rb = b[0]*rabs(0, 1)+b[1]*rabs(0, 0);
        if (std::abs(d[2]*r(1, 2) - d[1]*r(2, 2)) > ra+rb)
            return false;
    }
    {
        Real ra = a[0]*rabs(2, 0)+a[2]*rabs(0, 0);
        Real rb = b[1]*rabs(1, 2)+b[2]*rabs(1, 1);
        if (std::abs(d[0]*r(2, 0) - d[2]*r(0, 0)) > ra+rb)
            return false;
    }
    {
        Real ra = a[0]*rabs(2, 1)+a[2]*rabs(0, 1);
        Real rb = b[0]*rabs(1, 2)+b[2]*rabs(1, 0);
        if (std::abs(d[0]*r(2, 1) - d[2]*r(0, 1)) > ra+rb)
            return false;
    }
    {
        Real ra = a[0]*rabs(2, 2)+a[2]*rabs(0, 2);
        Real rb = b[0]*rabs(1, 1)+b[1]*rabs(1, 0);
        if (std::abs(d[0]*r(2, 2) - d[2]*r(0, 2)) > ra+rb)
            return false;
    }

    {
        Real ra = a[0]*rabs(1, 0)+a[1]*rabs(0, 0);
        Real rb = b[1]*rabs(2, 2)+b[2]*rabs(2, 1);
        if (std::abs(d[1]*r(0, 0) - d[0]*r(1, 0)) > ra+rb)
            return false;
    }
    {
        Real ra = a[0]*rabs(1, 1)+a[1]*rabs(0, 1);
        Real rb = b[0]*rabs(2, 2)+b[2]*rabs(2, 0);
        if (std::abs(d[1]*r(0, 1) - d[0]*r(1, 1)) > ra+rb)
            return false;
    }
    {
        Real ra = a[0]*rabs(1, 2)+a[1]*rabs(0, 2);
        Real rb = b[0]*rabs(2, 1)+b[1]*rabs(2, 0);
        if (std::abs(d[1]*r(0, 2) - d[0]*r(1, 2)) > ra+rb)
            return false;
    }
    return true;
}

bool OrientedBoundingBox::intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance) const {
    // Transform the ray to the bounding box's reference frame.
    
    Vec3 orig = ~getTransform()*origin;
    UnitVec3 dir = ~getTransform().R()*direction;
    
    // Check it against each plane that defines a side of the box.
    
    Real minDist = MostNegativeReal;
    Real maxDist = MostPositiveReal;
    if (dir[0] == 0.0) {
        if (orig[0] < 0 || orig[0] > getSize()[0])
            return false;
    }
    else {
        Real dist1 = -orig[0]/dir[0];
        Real dist2 = (getSize()[0]-orig[0])/dir[0];
        if (dist1 < dist2) {
            if (dist1 > minDist)
                minDist = dist1;
            if (dist2 < maxDist)
                maxDist = dist2;
      }
      else {
          if (dist2 > minDist)
              minDist = dist2;
          if (dist1 < maxDist)
              maxDist = dist1;
      }
      if (minDist > maxDist || maxDist < 0.0)
          return false;
    }
    if (dir[1] == 0.0) {
        if (orig[1] < 0 || orig[1] > getSize()[1])
            return false;
    }
    else {
        Real dist1 = -orig[1]/dir[1];
        Real dist2 = (getSize()[1]-orig[1])/dir[1];
        if (dist1 < dist2) {
            if (dist1 > minDist)
                minDist = dist1;
            if (dist2 < maxDist)
              maxDist = dist2;
        }
        else {
          if (dist2 > minDist)
              minDist = dist2;
          if (dist1 < maxDist)
              maxDist = dist1;
        }
        if (minDist > maxDist || maxDist < 0.0)
            return false;
    }
    if (dir[2] == 0.0) {
        if (orig[2] < 0 || orig[2] > getSize()[2])
            return false;
    }
    else {
        Real dist1 = -orig[2]/dir[2];
        Real dist2 = (getSize()[2]-orig[2])/dir[2];
        if (dist1 < dist2) {
            if (dist1 > minDist)
                minDist = dist1;
            if (dist2 < maxDist)
                maxDist = dist2;
        }
        else {
            if (dist2 > minDist)
                minDist = dist2;
            if (dist1 < maxDist)
                maxDist = dist1;
        }
        if (minDist > maxDist || maxDist < 0.0)
            return false;
    }
    if (minDist > 0)
        distance = minDist;
    else
        distance = 0;
    return true;
}

Vec3 OrientedBoundingBox::findNearestPoint(const Vec3& position) const {
    // Transform the point to the bounding box's reference frame.
    
    Vec3 p = ~getTransform()*position;
    
    // Find the nearest point in the box.
    
    if (p[0] < 0)
        p[0] = 0;
    if (p[0] > getSize()[0])
        p[0] = getSize()[0];
    if (p[1] < 0)
        p[1] = 0;
    if (p[1] > getSize()[1])
        p[1] = getSize()[1];
    if (p[2] < 0)
        p[2] = 0;
    if (p[2] > getSize()[2])
        p[2] = getSize()[2];
    
    // Transform it back again.
    
    return getTransform()*p;
}

void OrientedBoundingBox::getCorners(Vec3 corners[8]) const {
    Vec3 dx = size[0]*transform.R().col(0);
    Vec3 dy = size[1]*transform.R().col(1);
    Vec3 dz = size[2]*transform.R().col(2);
    corners[0] = transform.p();
    corners[1] = corners[0]+dx;
    corners[2] = corners[0]+dy;
    corners[3] = corners[1]+dy;
    corners[4] = corners[0]+dz;
    corners[5] = corners[1]+dz;
    corners[6] = corners[2]+dz;
    corners[7] = corners[3]+dz;
}

OrientedBoundingBox operator*(const Transform& t, const OrientedBoundingBox& box) {
    return OrientedBoundingBox(t*box.getTransform(), box.getSize());
}

} // namespace SimTK

