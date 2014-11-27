/* -------------------------------------------------------------------------- *
*                      Simbody(tm): Geometry Playground                      *
* -------------------------------------------------------------------------- *
* This is part of the SimTK biosimulation toolkit originating from           *
* Simbios, the NIH National Center for Physics-Based Simulation of           *
* Biological Structures at Stanford, funded under the NIH Roadmap for        *
* Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
*                                                                            *
* Portions copyright (c) 2007-12 Stanford University and the Authors.        *
* Authors: Michael Sherman                                                   *
* Contributors:                                                              *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License"); you may    *
* not use this file except in compliance with the License. You may obtain a  *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
* -------------------------------------------------------------------------- */

#include "Simbody.h"
#include "simmath/internal/OBBTree.h"

using namespace SimTK;
using std::cout; using std::endl;

// Trim r to 2 decimal places.
static Real trim2(Real r) {
	return std::floor(r * 100 + 0.5) / 100;
}

template <class P>
static void drawBoundingSphere(const Array_<Vec<3, P> >& v,
	const Transform& pose,
	MobilizedBody& body,
	const Vec3& color = Orange)
{
	Array_<int> which, support;
	Real tstart = realTime();
	Geo::Sphere_<P> pms = Geo::Point_<P>::calcBoundingSphere(v, which);
	printf("Time for %d points: %gs\n", v.size(), realTime() - tstart);
	cout << "ctr=" << pms.getCenter() << " rad=" << pms.getRadius() << "\n";
	cout << "volume=" << pms.findVolume() << "\n";
	for (unsigned i = 0; i<v.size(); ++i) {
		const P d = pms.getRadius() - (v[i] - pms.getCenter()).norm();
		if (d < 0)
			printf("*** %d outside by %.17g\n", i, d);
	}
	std::cout << "which=" << which << "\n";
	Decorations dec;
	//dec.addDecoration(Vec3(pms.getCenter()), 
	//    DecorativeSphere(pms.getRadius()).setColor(color).setResolution(5));
	//for (unsigned i=0; i < which.size(); ++i)
	//    dec.addDecoration(Vec3(v[which[i]]), DecorativePoint().setColor(Green).setScale(3));
	//body.addBodyDecoration(pose,
	//    dec.setOpacity(.3).setLineThickness(1));

	tstart = realTime();
	SimTK::Geo::OrientedBox_<P> obb = Geo::Point_<P>::calcOrientedBoundingBox(v, support);
	printf("OBB time for %d points: %gs\n", v.size(), realTime() - tstart);
	cout << "ctr=" << obb.getCenter() << " rad=" << obb.getHalfLengths() << "\n";
	cout << "volume=" << obb.getBox().findVolume() << "\n";


	Decorations obbDec;
	Transform X_FB; X_FB.updP() = Vec3(obb.getTransform().p());
	UnitVec3 x(Vec3(obb.getTransform().x())), y(Vec3(obb.getTransform().y()));
	X_FB.updR() = Rotation(x, XAxis, y, YAxis);
	obbDec.addDecoration(X_FB,
		DecorativeBrick(Vec3(obb.getHalfLengths())).setColor(Blue)
		//.setRepresentation(DecorativeGeometry::DrawWireframe)
		.setLineThickness(2).setOpacity(0.1));
	for (unsigned i = 0; i < support.size(); ++i) {
		const Vec<3, P>& supP = v[support[i]];
		const Vec3 sup((Real)supP[0], (Real)supP[1], (Real)supP[2]);
		obbDec.addDecoration(sup,
			DecorativePoint().setColor(Blue).setScale(4));
	}

	Array_<int> badSupport;
	Geo::OrientedBox_<P> badObb =
		Geo::Point_<P>::calcOrientedBoundingBox(v, badSupport, false);
	X_FB.updP() = Vec3(badObb.getTransform().p());
	x = UnitVec3(Vec3(badObb.getTransform().x()));
	y = UnitVec3(Vec3(badObb.getTransform().y()));
	X_FB.updR() = Rotation(x, XAxis, y, YAxis);
	obbDec.addDecoration(X_FB,
		DecorativeBrick(Vec3(badObb.getHalfLengths())).setColor(Gray)
		//.setRepresentation(DecorativeGeometry::DrawWireframe)
		.setOpacity(0.1));

	body.addBodyDecoration(pose, obbDec);
}


template <class P>
static void drawApproxBoundingSphere(const Array_<Vec<3, P> >& v,
	const Transform& pose,
	MobilizedBody& body,
	const Vec3& color = Gray)
{
	return;
	Real tstart = realTime();
	SimTK::Geo::Sphere_<P> pms = Geo::Point_<P>::calcApproxBoundingSphere(v);
	printf("Time for approx sphere around %d points: %gs\n",
		v.size(), realTime() - tstart);
	cout << "ctr=" << pms.getCenter() << " rad=" << pms.getRadius() << "\n";
	cout << "volume=" << pms.findVolume() << "\n";
	for (unsigned i = 0; i<v.size(); ++i) {
		const P d = pms.getRadius() - (v[i] - pms.getCenter()).norm();
		if (d < 0)
			printf("*** %d outside by %.17g\n", i, d);
	}
	Decorations dec;
	dec.addDecoration(Vec3(pms.getCenter()),
		DecorativeSphere(pms.getRadius()).setColor(color).setResolution(5));
	body.addBodyDecoration(pose,
		dec.setOpacity(.1).setLineThickness(1));
}

static void draw(const Geo::CubicBezierCurve& curve, const Transform& pose,
	MobilizedBody& body, const Vec3& color = Red, bool control = true)
{
	const int Resolution = 10;
	Decorations dec; dec.setColor(color).setLineThickness(2);

	Vec3 prev = curve.evalP(0);
	for (int i = 1; i <= Resolution; ++i) {
		const Real u = (Real)i / Resolution;
		Vec3 p = curve.evalP(u);
		dec.addDecoration(DecorativeLine(prev, p));
		prev = p;
	}
	if (control) {
		Vec<4, Vec3> B = curve.getControlPoints();
		for (int i = 0; i<4; ++i) {
			dec.addDecoration(DecorativePoint(B[i]));
			dec.addDecoration(DecorativeLine(B[i], B[(i + 1) % 4])
				.setColor(Gray).setLineThickness(1));
		}
		dec.addDecoration(DecorativeLine(B[0], B[2])
			.setColor(Green).setLineThickness(1));
		dec.addDecoration(DecorativeLine(B[1], B[3])
			.setColor(Green).setLineThickness(1));

		const int NumFrames = 0 * 100;
		for (int i = 0; i < NumFrames; ++i) {
			const Real u = (Real)i / (NumFrames - 1);
			Transform X_FP;
			const Real k = curve.calcCurveFrame(u, X_FP);
			dec.addDecoration(
				DecorativeLine(X_FP.p(), X_FP.p() + X_FP.x() / 4)
				.setColor(Blue).setLineThickness(1));
			dec.addDecoration(
				DecorativeLine(X_FP.p(), X_FP.p() + X_FP.z() / 8)
				.setColor(Red).setLineThickness(1));
		}
	}
	body.addBodyDecoration(pose, dec);
}

static void drawLines(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
	const Transform& pose,
	MobilizedBody& body, const Vec3& color = Black)
{
	body.addBodyDecoration(pose, DecorativeLine(a, b).setLineThickness(2));
	body.addBodyDecoration(pose, DecorativeLine(b, c).setLineThickness(2));
	body.addBodyDecoration(pose, DecorativeLine(c, d).setLineThickness(2));
	body.addBodyDecoration(pose, DecorativeLine(d, a).setLineThickness(2));
}


static void draw(const Geo::BicubicBezierPatch& patch, const Transform& pose,
	MobilizedBody& body, const Vec3& color = Red, bool control = true)
{
	const int Resolution = 10;
	Mat<4, 4, Vec3> B = patch.getControlPoints();

	for (int i = 0; i <= Resolution; ++i) {
		const Real t = Real(i) / Resolution;
		Geo::CubicBezierCurve iso = Geo::BicubicBezierPatch::calcIsoCurveU(B, t);
		draw(iso, pose, body, color, false);
		iso = Geo::BicubicBezierPatch::calcIsoCurveW(B, t);
		draw(iso, pose, body, color, false);
	}

	if (control) {
		for (int i = 0; i<4; ++i)
		for (int j = 0; j<4; ++j) {
			body.addBodyDecoration(pose,
				DecorativePoint(B[i][j]).setScale(4).setColor(Magenta));
		}
		drawLines(B[0][0], B[0][3], B[3][3], B[3][0], pose, body);//base
		drawLines(B[1][1], B[1][2], B[2][2], B[2][1], pose, body);//cap

		drawLines(B[0][0], B[0][1], B[0][2], B[0][3], pose, body);//bounds
		drawLines(B[3][0], B[3][1], B[3][2], B[3][3], pose, body);
		drawLines(B[0][0], B[1][0], B[2][0], B[3][0], pose, body);
		drawLines(B[0][3], B[1][3], B[2][3], B[3][3], pose, body);

		drawLines(B[0][0], B[0][1], B[1][1], B[1][0], pose, body);//corners
		drawLines(B[2][0], B[2][1], B[3][1], B[3][0], pose, body);
		drawLines(B[2][2], B[2][3], B[3][3], B[3][2], pose, body);
		drawLines(B[0][2], B[0][3], B[1][3], B[1][2], pose, body);

		drawLines(B[0][1], B[0][2], B[1][2], B[1][1], pose, body);//middles
		drawLines(B[1][0], B[1][1], B[2][1], B[2][0], pose, body);
		drawLines(B[2][1], B[2][2], B[3][2], B[3][1], pose, body);
		drawLines(B[1][2], B[1][3], B[2][3], B[2][2], pose, body);

		Geo::OrientedBox bs = patch.calcOrientedBoundingBox();
		body.addBodyDecoration(pose*bs.getTransform(),
			DecorativeBrick(bs.getHalfLengths()).setOpacity(.2)
			.setResolution(10));
	}
}

void makeDecorations(int level, const OBBNode& node, int which,
	Decorations& dec) {
	static Vec3 colors[6] = { Cyan, Green, Magenta, Blue, Orange, Red };

	Vec3 color = colors[node.height % 6];

	if (node.height == 0 || node.depth<2)
		dec.addDecoration(node.box.getTransform(),
		DecorativeBrick(node.box.getHalfLengths()).setColor(color)
		.setOpacity(Real(.9) / (node.height + 1)));

	for (unsigned i = 0; i<node.children.size(); ++i)
		makeDecorations(level + 1, node.children[i], which, dec);
}

void draw(const OBBTree& tree, int which,
	const Transform& pose, MobilizedBody& body) {
	const OBBNode& root = tree.getRoot();
	Decorations dec;
	makeDecorations(0, root, which, dec);
	body.addBodyDecoration(pose, dec);
}

void drawSpline(const Vector& x, const Vector& y,
	MobilizedBody& body) {

	Vec3 colors[3] = { Blue, Magenta, Green };

	for (int m = 1; m <= 3; ++m) {
		int order = 2 * m - 1;

		Spline_<Real> curve[3] = {
			SplineFitter<Real>::fitForSmoothingParameter(order, x, y, 0)
			.getSpline(),
			SplineFitter<Real>::fitForSmoothingParameter(order, x, y, .1)
			.getSpline(),
			SplineFitter<Real>::fitForSmoothingParameter(order, x, y, 10)
			.getSpline()
		};

		Vec3 offs[3] = { Vec3(10 * m, 0, 0), Vec3(10 * m, -3, 0), Vec3(10 * m, -6, 0) };
		for (int c = 0; c<3; ++c)
		for (int i = 0; i < x.size(); ++i) {
			body.addBodyDecoration(offs[c],
				DecorativePoint(Vec3(x[i], y[i], 0)).setScale(2)
				.setColor(colors[c]));
		}

		const int NSegs = 100;
		const Real dx = (x[x.size() - 1] - x[0]) / NSegs;
		for (int c = 0; c<3; ++c)
		for (int s = 0; s<NSegs; ++s) {
			Real xx = x[0] + s*dx;
			Vec3 p1(xx, curve[c].calcValue(Vector(1, xx)), 0);
			Vec3 p2(xx + dx, curve[c].calcValue(Vector(1, xx + dx)), 0);
			body.addBodyDecoration(offs[c],
				DecorativeLine(p1, p2).setColor(colors[c])
				.setLineThickness(3));
		}

	}
}



int main() {

	// Create the system, with subsystems for the bodies and some forces.
	MultibodySystem system;
	system.setUseUniformBackground(true);
	SimbodyMatterSubsystem matter(system);
	GeneralForceSubsystem forces(system);

	Real spxdata[] = { 0, 1, 2, 3, 4, 5, 6 };
	Real spydata[] = { 0, 1, 1.5, 2, 1, 1, 3 };
	Vector spx(7, spxdata), spy(7, spydata);
	drawSpline(spx, spy, matter.Ground());


	// SuperEllipsoids 
	// -------------------------------------------------------------------------------
	ContactGeometry::SuperEllipsoid SE_geom_01(Vec3(0.5, 0.5, 0.5) /*dimensions*/, Vec2(2.0,2.0) /*exponents*/);
	matter.updGround().addBodyDecoration(Transform(Vec3(15.0, 40.0, 0.0)), SE_geom_01.createDecorativeGeometry()
		.setColor(Red)
		.setOpacity(0.9));
	matter.Ground().addBodyDecoration(Transform(Vec3(15.0, 40.0, 0.0)), SE_geom_01.createDecorativeGeometry()
		.setRepresentation(DecorativeGeometry::DrawWireframe)
		.setColor(Black));

	ContactGeometry::SuperEllipsoid SE_geom_02(Vec3(0.5, 0.5, 1.0), Vec2(3.0, 3.0));
	matter.updGround().addBodyDecoration(Transform(Vec3(16.5,40.0,0.0)), SE_geom_02.createDecorativeGeometry()
		.setColor(Green)
		.setOpacity(0.9));
	matter.Ground().addBodyDecoration(Transform(Vec3(16.5,40.0,0.0)), SE_geom_02.createDecorativeGeometry()
		.setRepresentation(DecorativeGeometry::DrawWireframe)
		.setColor(Black));

	ContactGeometry::SuperEllipsoid SE_geom_03(Vec3(0.5, 0.5, 1.5), Vec2(10.0, 10.0));
	matter.updGround().addBodyDecoration(Transform(Vec3(18.0,40.0,0.0)), SE_geom_03.createDecorativeGeometry()
		.setColor(Blue)
		.setOpacity(0.9));
	matter.Ground().addBodyDecoration(Transform(Vec3(18.0,40.0,0.0)), SE_geom_03.createDecorativeGeometry()
		.setRepresentation(DecorativeGeometry::DrawWireframe)
		.setColor(Black));

	// -------------------------------------------------------------------------------



	Vec<4, Vec3> B(Vec3(3, 0, 0), Vec3(5, .5, 2), Vec3(4, 1, 0), Vec3(6, 0, 0));
	Geo::CubicBezierCurve curve(B);
	draw(curve, Vec3(0), matter.Ground(), Orange);
	Geo::CubicBezierCurve left, right, rl, rr;
	curve.split(.25, left, right);
	right.bisect(rl, rr);
	Vec3 offs(0, 0, -2);
	draw(left, offs, matter.Ground(), Blue);
	draw(rl, offs, matter.Ground(), Green);
	draw(rr, offs, matter.Ground(), Magenta);
	Geo::OrientedBox leftOBB = left.calcOrientedBoundingBox();
	Geo::OrientedBox rlOBB = rl.calcOrientedBoundingBox();
	Geo::OrientedBox rrOBB = rr.calcOrientedBoundingBox();
	matter.Ground().addBodyDecoration(leftOBB.getTransform() + offs,
		DecorativeBrick(leftOBB.getHalfLengths()).setOpacity(.1).setColor(Blue));
	matter.Ground().addBodyDecoration(rlOBB.getTransform() + offs,
		DecorativeBrick(rlOBB.getHalfLengths()).setOpacity(.1).setColor(Green));
	matter.Ground().addBodyDecoration(rrOBB.getTransform() + offs,
		DecorativeBrick(rrOBB.getHalfLengths()).setOpacity(.1).setColor(Magenta));

	Geo::Box box(Vec3(3, 4, 2)); // half lengths

	Geo::OrientedBox obox(Transform(), Vec3(1, 2, 3));
	// Non-intersecting box for which no face will serve as separator.
	// In this case the fast method can't tell they are separated.
	obox.setTransform(Transform(
		Rotation(BodyRotationSequence, Pi / 4, XAxis, Pi / 8, YAxis, -Pi / 4, ZAxis),
		Vec3(1.5, -5, 5.25)));
	matter.Ground().addBodyDecoration(Vec3(0),
		DecorativeBrick(box.getHalfLengths()).setOpacity(.1).setColor(Blue));
	matter.Ground().addBodyDecoration(obox.getTransform(),
		DecorativeBrick(obox.getHalfLengths()).setOpacity(.1).setColor(Green));
	cout << "May intersect=" << box.mayIntersectOrientedBox(obox) << "\n";
	cout << "Intersects=" << box.intersectsOrientedBox(obox) << "\n";

	Geo::Sphere curveSphere = curve.calcBoundingSphere();
	Geo::AlignedBox curveAABB = curve.calcAxisAlignedBoundingBox();
	Geo::OrientedBox curveOBB = curve.calcOrientedBoundingBox();
	cout << "vs=" << curveSphere.findVolume()
		<< " as=" << curveAABB.getBox().findVolume()
		<< " os=" << curveOBB.getBox().findVolume() << "\n";
	//matter.Ground().addBodyDecoration(curveSphere.getCenter(),
	//    DecorativeSphere(curveSphere.getRadius()).setOpacity(.1));
	//matter.Ground().addBodyDecoration(curveAABB.getCenter(),
	//    DecorativeBrick(curveAABB.getHalfLengths())
	//        .setRepresentation(DecorativeGeometry::DrawWireframe)
	//        .setColor(Gray));
	//matter.Ground().addBodyDecoration(curveOBB.getTransform(),
	//    DecorativeBrick(curveOBB.getHalfLengths())
	//        .setOpacity(.1)
	//        .setColor(Blue));
	const int n = 5;
	const fVec3 shft(0, 0, 0);
	fVec3 pd[] = { shft + fVec3(.5f, 0, -.5f), shft + fVec3(0, 1, .4f), shft + fVec3(-.5f, 1e-3f, 0)
		, shft + fVec3(.5f, -.2f, .1f), shft + fVec3(.5f, .5f, .5f)
	};
	Array_<int> which;
	//Array_<fVec3> p(pd, pd+n);
	Array_<fVec3> p;
	fVec3 start(0, 0, 0); float r = .7f; float noise = 1e-7f;
	p.push_back(fVec3(r, 0, 0));
	p.push_back(fVec3(0, r, 0));
	p.push_back(fVec3(0, 0, r));
	for (int i = 0; i<n - 3; ++i) {
		UnitVec3 uv(Test::randVec3());
		float nz = noise*(float)Test::randReal();
		fVec3 fuv((float)uv[0], (float)uv[1], (float)uv[2]);
		fVec3 fru = start + (r + nz)*fuv;
		p.push_back(fru);
	}
	drawBoundingSphere<float>(p, Transform(), matter.Ground());
	drawApproxBoundingSphere<float>(p, Transform(), matter.Ground());

	const int Nx = 4, Ny = 5;
	const Real xData[Nx] = { .1, 1, 2, 4 };
	const Real yData[Ny] = { -3, -2, 0, 1, 3 };
	const Real fData[Nx*Ny] = { 1, 2, 3, 3, 2,
		1.1, 2.1, 3.1, 3.1, 2.1,
		1, 2, 7, 3, 2,
		1.2, 2.2, 3.2, 3.2, 2.2 };
	const Vector x(Nx, xData);
	const Vector y(Ny, yData);
	const Matrix f(Nx, Ny, fData);
	BicubicSurface rough(x, y, f, 0);
	BicubicSurface smooth(x, y, f, 1);



	Vector xp(Vec2(.25, 3.25)), yp(Vec2(.75, 5.75));
	// Worst possible patch?:
	Matrix fp(Mat22(1.5, 1.7,
		1.3, 1.6));
	Matrix fxp(2 * Mat22(-1.1, -1.2,
		-1.3, -1.4));
	Matrix fyp(-2 * Mat22(1.2, 1.3,
		1.4, 1.1));
	Matrix fxyp(2 * Mat22(1.01, -1.02,
		-1.04, 1.03));
	// Nice patch:
	//Matrix fp(Mat22(1, 1,
	//                1, 1));
	//Matrix fxp(.5*Mat22(-1,  0,
	//                    0, -1));
	//Matrix fyp(.5*Mat22(-1,  0,
	//                    0,  -1));
	//Matrix fxyp(0*Mat22(1, -1,
	//                    -1, 1));
	// One-hump patch:
	//Matrix fp(Mat22(1, 1,
	//                1, 1));
	//Matrix fxp(Mat22(1,  1,
	//                 -1, -1));
	//Matrix fyp(Mat22(1,  -1,
	//                 1,  -1));
	//Matrix fxyp(0.5*Mat22(1, 3,
	//                    -3, 4));

	BicubicSurface patch(xp, yp, fp, fxp, fyp, fxyp);
	Rotation xm90(-Pi / 2, XAxis);
	Transform patchPose(xm90, Vec3(4, 2, 0));


	int nx, ny; patch.getNumPatches(nx, ny);
	printf("surface 'patch' has %dx%d patches\n", nx, ny);
	Geo::BicubicBezierPatch bpatch = patch.calcBezierPatch(0, 0);
	const Mat<4, 4, Vec3>& pts = bpatch.getControlPoints();
	cout << "Bezier pts:\n" << pts;
	Geo::BicubicBezierPatch patch00, patch01, patch10, patch11;
	bpatch.split(.2, .3, patch00, patch01, patch10, patch11);

	draw(patch00, patchPose + Vec3(5, 0, 0), matter.Ground(), Cyan);
	draw(patch01, patchPose + Vec3(5, 0, 0), matter.Ground(), Green);
	draw(patch10, patchPose + Vec3(5, 0, 0), matter.Ground(), Purple);
	draw(patch11, patchPose + Vec3(5, 0, 0), matter.Ground(), Blue);



	for (int i = 0; i<4; ++i)
	for (int j = 0; j<4; ++j)
		matter.Ground().addBodyDecoration(patchPose,
		DecorativePoint(pts[i][j]).setScale(4).setColor(Magenta));

	draw(bpatch.getBoundaryCurveU0(), patchPose, matter.Ground());
	draw(bpatch.calcIsoCurveU(pts, .5), patchPose, matter.Ground());
	draw(bpatch.getBoundaryCurveU1(), patchPose, matter.Ground());
	draw(bpatch.getBoundaryCurveW0(), patchPose, matter.Ground(), Blue);
	draw(bpatch.calcIsoCurveW(pts, .5), patchPose, matter.Ground(), Blue);
	draw(bpatch.getBoundaryCurveW1(), patchPose, matter.Ground(), Blue);

	Real resolution = 31;
	PolygonalMesh patchMesh = patch.createPolygonalMesh(resolution);

	Array_<Vec3> v(patchMesh.getNumVertices());
	for (int i = 0; i < patchMesh.getNumVertices(); ++i)
		v[i] = patchMesh.getVertexPosition(i);
	//drawBoundingSphere<Real>(v, patchPose, matter.Ground());
	//drawApproxBoundingSphere<Real>(v, patchPose, matter.Ground());

	matter.Ground().addBodyDecoration(patchPose,
		DecorativeMesh(patchMesh).setRepresentation(DecorativeGeometry::DrawWireframe)
		.setColor(Gray));
	//matter.Ground().addBodyDecoration(patchPose, 
	//    DecorativeMesh(patchMesh)
	//    .setColor(Blue).setOpacity(.4));

	matter.Ground().addBodyDecoration(patchPose,
		Decorations()
		.addDecoration(DecorativePoint(Vec3(xp[0], yp[0], fp(0, 0))))
		.addDecoration(DecorativePoint(Vec3(xp[0], yp[1], fp(0, 1))).setColor(Blue))
		.addDecoration(DecorativePoint(Vec3(xp[1], yp[0], fp(1, 0))))
		.addDecoration(DecorativePoint(Vec3(xp[1], yp[1], fp(1, 1))))
		);

	//BicubicSurface surf(Vec2(0,0), Vec2(1,1), f, 0);
	Array_<int> dxx(2, 0), dyy(2, 0), dxy(2);
	dxy[0] = 0; dxy[1] = 1;
	Array_<int> dx(1, 0), dy(1, 1);

	Transform pose1(xm90, Vec3(0, -15, 0));
	Transform pose2(xm90, Vec3(5, -15, 0));


	ContactGeometry::SmoothHeightMap smoothMap(smooth);

	draw(smoothMap.getOBBTree(), 1, pose2, matter.Ground());

	//#ifdef NOTDEF
	for (int s = 0; s <= 1; ++s) {
		BicubicSurface surf = s ? smooth : rough;
		Transform pose = s ? pose2 : pose1;
		for (int i = 0; i<Nx; ++i)
		for (int j = 0; j<Ny; ++j) {
			const Vec2 pt(x[i], y[j]);

			Vec2 k; Transform X_SP;
			surf.calcParaboloid(pt, X_SP, k);
			const Vec3& P = X_SP.p();

			// The original tangents.
			const Vec3 dPdx(1, 0, surf.calcDerivative(dx, pt));
			const Vec3 dPdy(0, 1, surf.calcDerivative(dy, pt));

			const Real BigLen = 0.5, TextScale = .1;
			matter.Ground().addBodyDecoration(pose,
				DecorativeLine(P, P + BigLen*X_SP.z()).setColor(Black));
			matter.Ground().addBodyDecoration(pose,
				DecorativeLine(P, P + BigLen*X_SP.x()).setColor(Red).setLineThickness(3));
			matter.Ground().addBodyDecoration(
				Transform(pose.p() + pose.R()*(P + 1.1*BigLen*X_SP.x())),
				DecorativeText(String(trim2(k[0]))).setColor(Red).setScale(TextScale)
				.setFaceCamera(false));
			matter.Ground().addBodyDecoration(pose,
				DecorativeLine(P, P + BigLen*X_SP.y()).setColor(Blue).setLineThickness(3));
			matter.Ground().addBodyDecoration(
				Transform(pose.p() + pose.R()*(P + 1.1*BigLen*X_SP.y())),
				DecorativeText(String(trim2(k[1]))).setColor(Blue).setScale(TextScale));
			matter.Ground().addBodyDecoration(pose,
				DecorativeLine(P, P + .75*BigLen*UnitVec3(dPdx)).setColor(Red));
			matter.Ground().addBodyDecoration(pose,
				DecorativeLine(P, P + .75*BigLen*UnitVec3(dPdy)).setColor(Blue));
		}
	}

	PolygonalMesh origMesh = rough.createPolygonalMesh(0);
	PolygonalMesh roughMesh = rough.createPolygonalMesh(resolution);
	PolygonalMesh smoothMesh = smooth.createPolygonalMesh(resolution);
	printf("stats: access=%d, same pt=%d, same patch=%d, nearby=%d\n",
		rough.getNumAccesses(), rough.getNumAccessesSamePoint(),
		rough.getNumAccessesSamePatch(), rough.getNumAccessesNearbyPatch());

	v.resize(roughMesh.getNumVertices());
	for (int i = 0; i < roughMesh.getNumVertices(); ++i)
		v[i] = roughMesh.getVertexPosition(i);
	drawBoundingSphere(v, pose1, matter.Ground());
	drawApproxBoundingSphere(v, pose1, matter.Ground());

	v.resize(smoothMesh.getNumVertices());
	for (int i = 0; i < smoothMesh.getNumVertices(); ++i)
		v[i] = smoothMesh.getVertexPosition(i);
	//drawBoundingSphere(v, pose2, matter.Ground());
	//drawApproxBoundingSphere(v, pose2, matter.Ground());

	v.resize(smoothMesh.getNumVertices() + roughMesh.getNumVertices());
	const Transform X_21 = ~pose2*pose1;
	for (int i = 0; i < roughMesh.getNumVertices(); ++i)
		v[smoothMesh.getNumVertices() + i] = X_21*roughMesh.getVertexPosition(i);
	//drawBoundingSphere(v, pose2, matter.Ground(), Blue);
	//drawApproxBoundingSphere(v, pose2, matter.Ground(), Blue);

	//matter.Ground().addBodyDecoration(pose1, 
	//    DecorativeMesh(origMesh).setRepresentation(DecorativeGeometry::DrawWireframe)
	//    .setColor(Gray));
	matter.Ground().addBodyDecoration(pose1,
		DecorativeMesh(roughMesh)
		.setColor(Gray).setOpacity(.4));

	//matter.Ground().addBodyDecoration(pose2, 
	//    DecorativeMesh(origMesh).setRepresentation(DecorativeGeometry::DrawWireframe)
	//    .setColor(Gray));
	matter.Ground().addBodyDecoration(pose2,
		DecorativeMesh(smoothMesh)
		.setColor(Black).setOpacity(.4)
		.setRepresentation(DecorativeGeometry::DrawWireframe));

	BicubicSurface::PatchHint hint1, hint2;
	Array_<int> fx(1, 0), fy(1, 1);
	const Real Len = .25;
	Vec2 xy0 = rough.getMinXY(), range = rough.getMaxXY() - xy0;
	const int NSamples = 25;
	const Vec2 incr(range / (NSamples - 1));
	//for (int i=0; i<NSamples; ++i)
	//    for (int j=0; j<NSamples; ++j) {
	//        Vec2 xy(xy0 + Vec2(i*incr[0], j*incr[1]));
	//        Vec3 start1(xy[0],xy[1], rough.calcValue(xy,hint1));
	//        UnitVec3 nz1 = rough.calcUnitNormal(xy,hint1);
	//        matter.Ground().addBodyDecoration(pose1,
	//            DecorativeLine(start1, start1+Len*nz1)
	//                .setColor(Green));

	//        Vec3 start2(xy[0],xy[1], smooth.calcValue(xy,hint2));
	//        UnitVec3 nz2 = smooth.calcUnitNormal(xy,hint2);
	//        matter.Ground().addBodyDecoration(pose2,
	//            DecorativeLine(start2, start2+Len*nz2)
	//                .setColor(Green));
	//    }
	//#endif

	// Add gravity as a force element.
	Rotation x45(Pi / 4, XAxis);
	Rotation y45(Pi / 4, YAxis);
	Rotation z45(Pi / 4, ZAxis);
	Force::UniformGravity gravity(forces, matter, Vec3(10, -9.8, 3));
	// Create the body and some artwork for it.
	Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
	pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1).setColor(Red));
	pendulumBody.addDecoration(Transform(), DecorativePoint(Vec3(0, .5, 0)).setColor(Green));
	// Add an instance of the body to the multibody system by connecting
	// it to Ground via a pin mobilizer.
	MobilizedBody::Pin pendulum1(matter.updGround(),
		Transform(/*x45,*/Vec3(0, -1, 0)),
		pendulumBody,
		Transform(Vec3(0, 1, 0)));
	MobilizedBody::Pin pendulum1b(pendulum1,
		Transform(/*x45,*/Vec3(0, -1, 0)),
		pendulumBody,
		Transform(Vec3(0, 1, 0)));


	MobilizedBody::Free pendulum2(matter.updGround(),
		Transform(/*x45,*/Vec3(2, -1, 0)),
		pendulumBody,
		Transform(Vec3(0, 1, 0)));
	Constraint::Ball ballcons2(matter.updGround(), Vec3(2, -1, 0),
		pendulum2, Vec3(0, 1, 0));
	const Transform& X_GF2 = pendulum2.getDefaultInboardFrame();
	const Transform& X_P2M = pendulum2.getDefaultOutboardFrame();
	Constraint::ConstantAngle angx2(matter.Ground(), X_GF2.x(),
		pendulum2, X_P2M.z());
	Constraint::ConstantAngle angy2(matter.Ground(), X_GF2.y(),
		pendulum2, X_P2M.z());

	MobilizedBody::Free pendulum2b(pendulum2,
		Transform(/*x45,*/Vec3(0, -1, 0)),
		pendulumBody,
		Transform(Vec3(0, 1, 0)));
	Constraint::Ball ballcons2b(pendulum2, Vec3(0, -1, 0),
		pendulum2b, Vec3(0, 1, 0));
	const Transform& X_GF2b = pendulum2b.getDefaultInboardFrame();
	const Transform& X_P2Mb = pendulum2b.getDefaultOutboardFrame();
	Constraint::ConstantAngle angx2b(pendulum2, X_GF2b.x(),
		pendulum2b, X_P2Mb.z());
	Constraint::ConstantAngle angy2b(pendulum2, X_GF2b.y(),
		pendulum2b, X_P2Mb.z());


	// Visualize with default options; ask for a report every 1/30 of a second
	// to match the Visualizer's default 30 frames per second rate.
	Visualizer viz(system);
	system.addEventReporter(new Visualizer::Reporter(viz, 1. / 30));

	// Initialize the system and state.

	system.realizeTopology();
	State state = system.getDefaultState();
	pendulum1.setOneQ(state, 0, Pi / 4);
	//pendulum1.setOneU(state, 0, 1.0); // initial velocity 1 rad/sec

	//pendulum1b.setOneU(state, 0, 1.0); // initial velocity 1 rad/sec
	pendulum1b.setOneQ(state, 0, Pi / 4);

	pendulum2.setQToFitRotation(state, Rotation(Pi / 4, ZAxis));
	//pendulum2.setUToFitAngularVelocity(state, Vec3(0,0,1));
	pendulum2b.setQToFitRotation(state, Rotation(Pi / 4, ZAxis));
	//pendulum2b.setUToFitAngularVelocity(state, Vec3(0,0,1));

	system.realize(state);
	const Vector lambda = state.getMultipliers();
	Vector_<SpatialVec> consBodyForcesInG;
	Vector              consMobForces;
	matter.calcConstraintForcesFromMultipliers(state, -lambda, consBodyForcesInG,
		consMobForces);
	const MobodIndex p2x = pendulum2.getMobilizedBodyIndex();
	const MobodIndex p2bx = pendulum2b.getMobilizedBodyIndex();
	const Rotation& R_G2 = pendulum2.getBodyTransform(state).R();
	//consBodyForcesInG[p2x] = shiftForceFromTo(consBodyForcesInG[p2x],
	//                                         Vec3(0), R_G2*Vec3(0,1,0));

	const int nb = matter.getNumBodies();
	Vector_<SpatialVec> forcesAtMInG;
	matter.calcMobilizerReactionForces(state, forcesAtMInG);

	Vector_<SpatialVec> forcesAtFInG(nb); // to hold the result
	forcesAtFInG[0] = -forcesAtMInG[0]; // Ground is "welded" at origin
	for (MobilizedBodyIndex i(1); i < matter.getNumBodies(); ++i) {
		const MobilizedBody& body = matter.getMobilizedBody(i);
		const MobilizedBody& parent = body.getParentMobilizedBody();
		// Want to shift negated reaction by p_MF_G, the vector from M
		// to F across the mobilizer, expressed in Ground. We can get p_FM, 
		// then re-express in Ground for the shift and negate.
		const Vec3& p_FM = body.getMobilizerTransform(state).p();
		const Rotation& R_PF = body.getInboardFrame(state).R(); // In parent.
		const Rotation& R_GP = parent.getBodyTransform(state).R();
		Rotation R_GF = R_GP*R_PF;  // F frame orientation in Ground.
		Vec3     p_MF_G = -(R_GF*p_FM); // Re-express and negate shift vector. 
		forcesAtFInG[i] = -shiftForceBy(forcesAtMInG[i], p_MF_G);
	}

	std::cout << "Reactions @M: " << forcesAtMInG << "\n";
	std::cout << "Reactions @F: " << forcesAtFInG << "\n";

	const MobodIndex p1x = pendulum1.getMobilizedBodyIndex();
	const MobodIndex p1bx = pendulum1b.getMobilizedBodyIndex();
	const Rotation& R_G1 = pendulum1.getBodyTransform(state).R();
	const Rotation& R_G1b = pendulum1b.getBodyTransform(state).R();

	for (MobodIndex i(0); i < nb; ++i) {
		const Mobod& body = matter.getMobilizedBody(i);
		const Vec3&  p_BM = body.getOutboardFrame(state).p();
		const Rotation& R_GB = body.getBodyTransform(state).R();
		forcesAtMInG[i] = shiftForceFromTo(forcesAtMInG[i],
			R_GB*p_BM, Vec3(0));
	}


	std::cout << "FB_G=" << forcesAtMInG[p1x] << " " << forcesAtMInG[p1bx] << "\n";

	cout << "FC_G=" << -(ballcons2.getConstrainedBodyForcesAsVector(state)
		+ angx2.getConstrainedBodyForcesAsVector(state)
		+ angy2.getConstrainedBodyForcesAsVector(state))[1] << " ";
	cout << -(ballcons2b.getConstrainedBodyForcesAsVector(state)
		+ angx2b.getConstrainedBodyForcesAsVector(state)
		+ angy2b.getConstrainedBodyForcesAsVector(state))[1] << endl;

	viz.report(state);
	// Simulate it.
	getchar();

	RungeKuttaMersonIntegrator integ(system);
	TimeStepper ts(system, integ);
	ts.initialize(state);
	ts.stepTo(1.0);
	state = integ.getState();
	system.realize(state);
	matter.calcMobilizerReactionForces(state, forcesAtMInG);
	const Transform& X_GP = pendulum1.getBodyTransform(state);
	//forcesAtMInG[1][1] = X_GP.R()*forcesAtMInG[1][1];
	std::cout << "FM_G=" << forcesAtMInG << "\n";
	ts.stepTo(1.2);
	state = integ.getState();
	system.realize(state);
	matter.calcMobilizerReactionForces(state, forcesAtMInG);
	std::cout << "FM_G=" << forcesAtMInG << "\n";

}
