/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
 ****************************************************************************/
package jtoolex.voronoi;


import jtool.atom.IXYZ;
import jtool.atom.XYZ;

/**
 * Robust geometric predicates.
 * <p>
 * These geometric predicates are notoriously susceptible to roundoff error. For
 * example, the simplest and fastest test to determine whether a point c is left
 * of a line defined by two points a and b may fail when all three points are
 * nearly co-linear.
 * <p>
 * Therefore, each predicate is implemented by two types of methods. One method
 * is fast, but may yield incorrect answers. The other method is slower, because
 * it (1) computes a bound on the roundoff error and (2) reverts to an exact
 * algorithm if the fast method might yield the wrong answer.
 * <p>
 * Most applications should use the slower exact methods. The fast methods are
 * provided only for comparison.
 * <p>
 * These predicates are adapted from those developed by Jonathan Shewchuk, 1997,
 * Delaunay Refinement Mesh Generation: Ph.D. dissertation, Carnegie Mellon
 * University. (Currently, the methods here do not use Shewchuk's adaptive
 * four-stage pipeline. Instead, only two - the fastest and the exact stages -
 * are used.)
 * 
 * @author Dave Hale, Colorado School of Mines
 * @version 2001.04.03, 2006.08.02
 */
public class Geometry {
    
    /**
     * 为了可读性以及减少重复代码，这里使用会创建临时变量 {@link XYZ} 的
     * @author CHanzy
     * @return A, B, C 三点组成的三角形的面积，永远为正数
     */
    public static double area(IXYZ aA, IXYZ aB, IXYZ aC) {
        XYZ tA = XYZ.toXYZ(aA);
        XYZ tAB = aB.minus(tA);
        XYZ tAC = aC.minus(tA);
        tAB.cross2this(tAC);
        return 0.5 * tAB.norm();
    }
    /**
     * 确定点 D 是否位于由点 A、B 和 C 定义的平面的左侧。假定从平面的右侧看，ABC 满足逆时针的顺序。
     * @return 如果在左边则为正，右边则为负，刚好在平面上则为 0
     */
    public static double leftOfPlane(IXYZ aA, IXYZ aB, IXYZ aC, IXYZ aD) {
        return leftOfPlane(
            aA.x(), aA.y(), aA.z()
            , aB.x(), aB.y(), aB.z()
            , aC.x(), aC.y(), aC.z()
            , aD.x(), aD.y(), aD.z());
    }
    /**
     * 确定点 E 是否位于由点 A、B、C 和 D 定义的球体的内部。假定 {@code leftOfPlane(A, B, C, D) > 0}。
     * @return 如果在内部则为正，外部则为负，刚好在球面上则为 0
     */
    public static double inSphere(IXYZ aA, IXYZ aB, IXYZ aC, IXYZ aD, IXYZ aE) {
        return inSphere(
            aA.x(), aA.y(), aA.z()
            , aB.x(), aB.y(), aB.z()
            , aC.x(), aC.y(), aC.z()
            , aD.x(), aD.y(), aD.z()
            , aE.x(), aE.y(), aE.z());
    }
    /**
     * 计算由点 A，B，C 和 D 定义的球的中心。假定 {@code leftOfPlane(A, B, C, D) > 0}。
     * @return 球心 XYZ 坐标
     */
    public static XYZ centerSphere(IXYZ aA, IXYZ aB, IXYZ aC, IXYZ aD) {
        XYZ rCenter = new XYZ(0.0, 0.0, 0.0);
        centerSphere(
            aA.x(), aA.y(), aA.z()
            , aB.x(), aB.y(), aB.z()
            , aC.x(), aC.y(), aC.z()
            , aD.x(), aD.y(), aD.z()
            , rCenter);
        return rCenter;
    }
    
    
    /**
     * Computes the center of the sphere defined by the points a, b, c, and d. The
     * latter are assumed to be in CCW order, such that the method
     * {@link #leftOfPlane} would return a positive number.
     *
     * @param po array containing (x,y,z) coordinates of center.
     */
    public static void centerSphere(double xa, double ya, double za, double xb, double yb, double zb, double xc,
                                    double yc, double zc, double xd, double yd, double zd, XYZ po) {
        double adx = xa - xd;
        double bdx = xb - xd;
        double cdx = xc - xd;
        double ady = ya - yd;
        double bdy = yb - yd;
        double cdy = yc - yd;
        double adz = za - zd;
        double bdz = zb - zd;
        double cdz = zc - zd;
        double ads = adx * adx + ady * ady + adz * adz;
        double bds = bdx * bdx + bdy * bdy + bdz * bdz;
        double cds = cdx * cdx + cdy * cdy + cdz * cdz;
        double scale = 0.5 / leftOfPlane(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd);
        po.mX = xd
            + scale * (ads * (bdy * cdz - cdy * bdz) + bds * (cdy * adz - ady * cdz) + cds * (ady * bdz - bdy * adz));
        po.mY = yd
            + scale * (ads * (bdz * cdx - cdz * bdx) + bds * (cdz * adx - adz * cdx) + cds * (adz * bdx - bdz * adx));
        po.mZ = zd
            + scale * (ads * (bdx * cdy - cdx * bdy) + bds * (cdx * ady - adx * cdy) + cds * (adx * bdy - bdx * ady));
    }
    
    /**
     * Determines if a point e is inside the sphere defined by the points a, b, c,
     * and d. The latter are assumed to be in CCW order, such that the method
     * {@link #leftOfPlane} would return a positive number.
     *
     * @return positive, if inside the sphere; negative, if outside the sphere;
     *         zero, otherwise.
     */
    public static double inSphere(double xa, double ya, double za, double xb, double yb, double zb, double xc,
                                  double yc, double zc, double xd, double yd, double zd, double xe, double ye,
                                  double ze) {
        double aex = xa - xe;
        double bex = xb - xe;
        double cex = xc - xe;
        double dex = xd - xe;
        double aey = ya - ye;
        double bey = yb - ye;
        double cey = yc - ye;
        double dey = yd - ye;
        double aez = za - ze;
        double bez = zb - ze;
        double cez = zc - ze;
        double dez = zd - ze;
        
        double aexbey = aex * bey;
        double bexaey = bex * aey;
        double ab = aexbey - bexaey;
        double bexcey = bex * cey;
        double cexbey = cex * bey;
        double bc = bexcey - cexbey;
        double cexdey = cex * dey;
        double dexcey = dex * cey;
        double cd = cexdey - dexcey;
        double dexaey = dex * aey;
        double aexdey = aex * dey;
        double da = dexaey - aexdey;
        
        double aexcey = aex * cey;
        double cexaey = cex * aey;
        double ac = aexcey - cexaey;
        double bexdey = bex * dey;
        double dexbey = dex * bey;
        double bd = bexdey - dexbey;
        
        double abc = aez * bc - bez * ac + cez * ab;
        double bcd = bez * cd - cez * bd + dez * bc;
        double cda = cez * da + dez * ac + aez * cd;
        double dab = dez * ab + aez * bd + bez * da;
        
        double alift = aex * aex + aey * aey + aez * aez;
        double blift = bex * bex + bey * bey + bez * bez;
        double clift = cex * cex + cey * cey + cez * cez;
        double dlift = dex * dex + dey * dey + dez * dez;
        
        double det = dlift * abc - clift * dab + (blift * cda - alift * bcd);
        
        if (aez < 0.0) {
            aez = -aez;
        }
        if (bez < 0.0) {
            bez = -bez;
        }
        if (cez < 0.0) {
            cez = -cez;
        }
        if (dez < 0.0) {
            dez = -dez;
        }
        if (aexbey < 0.0) {
            aexbey = -aexbey;
        }
        if (bexaey < 0.0) {
            bexaey = -bexaey;
        }
        if (bexcey < 0.0) {
            bexcey = -bexcey;
        }
        if (cexbey < 0.0) {
            cexbey = -cexbey;
        }
        if (cexdey < 0.0) {
            cexdey = -cexdey;
        }
        if (dexcey < 0.0) {
            dexcey = -dexcey;
        }
        if (dexaey < 0.0) {
            dexaey = -dexaey;
        }
        if (aexdey < 0.0) {
            aexdey = -aexdey;
        }
        if (aexcey < 0.0) {
            aexcey = -aexcey;
        }
        if (cexaey < 0.0) {
            cexaey = -cexaey;
        }
        if (bexdey < 0.0) {
            bexdey = -bexdey;
        }
        if (dexbey < 0.0) {
            dexbey = -dexbey;
        }
        double permanent = ((cexdey + dexcey) * bez + (dexbey + bexdey) * cez + (bexcey + cexbey) * dez) * alift
            + ((dexaey + aexdey) * cez + (aexcey + cexaey) * dez + (cexdey + dexcey) * aez) * blift
            + ((aexbey + bexaey) * dez + (bexdey + dexbey) * aez + (dexaey + aexdey) * bez) * clift
            + ((bexcey + cexbey) * aez + (cexaey + aexcey) * bez + (aexbey + bexaey) * cez) * dlift;
        double errbound = INSERRBOUND * permanent;
        if (det > errbound || -det > errbound) {
            return det;
        }
        
        return inSphereExact(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd, xe, ye, ze);
    }
    
    
    
    /**
     * Determines if a point d is left of the plane defined by the points a, b, and
     * c. The latter are assumed to be in CCW order, as viewed from the right side
     * of the plane.
     *
     * @return positive, if left of plane; negative, if right of plane; zero,
     *         otherwise.
     */
    public static double leftOfPlane(double xa, double ya, double za, double xb, double yb, double zb, double xc,
                                     double yc, double zc, double xd, double yd, double zd) {
        double adx = xa - xd;
        double bdx = xb - xd;
        double cdx = xc - xd;
        double ady = ya - yd;
        double bdy = yb - yd;
        double cdy = yc - yd;
        double adz = za - zd;
        double bdz = zb - zd;
        double cdz = zc - zd;
        
        double bdxcdy = bdx * cdy;
        double cdxbdy = cdx * bdy;
        
        double cdxady = cdx * ady;
        double adxcdy = adx * cdy;
        
        double adxbdy = adx * bdy;
        double bdxady = bdx * ady;
        
        double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
        
        if (adz < 0.0) {
            adz = -adz;
        }
        if (bdz < 0.0) {
            bdz = -bdz;
        }
        if (cdz < 0.0) {
            cdz = -cdz;
        }
        if (bdxcdy < 0.0) {
            bdxcdy = -bdxcdy;
        }
        if (cdxbdy < 0.0) {
            cdxbdy = -cdxbdy;
        }
        if (cdxady < 0.0) {
            cdxady = -cdxady;
        }
        if (adxcdy < 0.0) {
            adxcdy = -adxcdy;
        }
        if (adxbdy < 0.0) {
            adxbdy = -adxbdy;
        }
        if (bdxady < 0.0) {
            bdxady = -bdxady;
        }
        double permanent = (bdxcdy + cdxbdy) * adz + (cdxady + adxcdy) * bdz + (adxbdy + bdxady) * cdz;
        double errbound = O3DERRBOUND * permanent;
        if (det > errbound || -det > errbound) {
            return det;
        }
        
        return leftOfPlaneExact(xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd);
    }
    
    
    private static class Two {
        double x, y;
    }
    private static final double EPSILON;
    private static final double INCERRBOUND;
    private static final double INSERRBOUND;
    private static final double IOSERRBOUND;
    private static final double O2DERRBOUND;
    private static final double O3DERRBOUND;
    private static final double SPLITTER;
    static {
        double epsilon = 1.0;
        double splitter = 1.0;
        boolean everyOther = true;
        do {
            epsilon *= 0.5;
            if (everyOther) {
                splitter *= 2.0;
            }
            everyOther = !everyOther;
        } while (1.0 + epsilon != 1.0);
        splitter += 1.0;
        EPSILON = epsilon;
        SPLITTER = splitter;
        O2DERRBOUND = 4.0 * EPSILON;
        O3DERRBOUND = 8.0 * EPSILON;
        INCERRBOUND = 11.0 * EPSILON;
        INSERRBOUND = 17.0 * EPSILON;
        IOSERRBOUND = 19.0 * EPSILON;
    }
    
    
    private static double inSphereExact(double xa, double ya, double za, double xb, double yb, double zb, double xc,
                                        double yc, double zc, double xd, double yd, double zd, double xe, double ye,
                                        double ze) {
        Two t = new Two();
        twoDiff(xa, xe, t);
        double aex = t.x;
        double aextail = t.y;
        twoDiff(ya, ye, t);
        double aey = t.x;
        double aeytail = t.y;
        twoDiff(za, ze, t);
        double aez = t.x;
        double aeztail = t.y;
        twoDiff(xb, xe, t);
        double bex = t.x;
        double bextail = t.y;
        twoDiff(yb, ye, t);
        double bey = t.x;
        double beytail = t.y;
        twoDiff(zb, ze, t);
        double bez = t.x;
        double beztail = t.y;
        twoDiff(xc, xe, t);
        double cex = t.x;
        double cextail = t.y;
        twoDiff(yc, ye, t);
        double cey = t.x;
        double ceytail = t.y;
        twoDiff(zc, ze, t);
        double cez = t.x;
        double ceztail = t.y;
        twoDiff(xd, xe, t);
        double dex = t.x;
        double dextail = t.y;
        twoDiff(yd, ye, t);
        double dey = t.x;
        double deytail = t.y;
        twoDiff(zd, ze, t);
        double dez = t.x;
        double deztail = t.y;
        
        double[] axby = new double[8];
        double[] bxay = new double[8];
        double[] ab = new double[16];
        twoTwoProduct(aex, aextail, bey, beytail, axby);
        double negate = -aey;
        double negatetail = -aeytail;
        twoTwoProduct(bex, bextail, negate, negatetail, bxay);
        int ablen = expansionSumZeroElimFast(8, axby, 8, bxay, ab);
        
        double[] bxcy = new double[8];
        double[] cxby = new double[8];
        double[] bc = new double[16];
        twoTwoProduct(bex, bextail, cey, ceytail, bxcy);
        negate = -bey;
        negatetail = -beytail;
        twoTwoProduct(cex, cextail, negate, negatetail, cxby);
        int bclen = expansionSumZeroElimFast(8, bxcy, 8, cxby, bc);
        
        double[] cxdy = new double[8];
        double[] dxcy = new double[8];
        double[] cd = new double[16];
        twoTwoProduct(cex, cextail, dey, deytail, cxdy);
        negate = -cey;
        negatetail = -ceytail;
        twoTwoProduct(dex, dextail, negate, negatetail, dxcy);
        int cdlen = expansionSumZeroElimFast(8, cxdy, 8, dxcy, cd);
        
        double[] dxay = new double[8];
        double[] axdy = new double[8];
        double[] da = new double[16];
        twoTwoProduct(dex, dextail, aey, aeytail, dxay);
        negate = -dey;
        negatetail = -deytail;
        twoTwoProduct(aex, aextail, negate, negatetail, axdy);
        int dalen = expansionSumZeroElimFast(8, dxay, 8, axdy, da);
        
        double[] axcy = new double[8];
        double[] cxay = new double[8];
        double[] ac = new double[16];
        twoTwoProduct(aex, aextail, cey, ceytail, axcy);
        negate = -aey;
        negatetail = -aeytail;
        twoTwoProduct(cex, cextail, negate, negatetail, cxay);
        int aclen = expansionSumZeroElimFast(8, axcy, 8, cxay, ac);
        
        double[] bxdy = new double[8];
        double[] dxby = new double[8];
        double[] bd = new double[16];
        twoTwoProduct(bex, bextail, dey, deytail, bxdy);
        negate = -bey;
        negatetail = -beytail;
        twoTwoProduct(dex, dextail, negate, negatetail, dxby);
        int bdlen = expansionSumZeroElimFast(8, bxdy, 8, dxby, bd);
        
        double[] t32a = new double[32];
        double[] t32b = new double[32];
        double[] t64a = new double[64];
        double[] t64b = new double[64];
        double[] t64c = new double[64];
        double[] t128 = new double[128];
        double[] t192 = new double[192];
        int t32alen, t32blen, t64alen, t64blen, t64clen, t128len, t192len;
        t32alen = scaleExpansionZeroElim(cdlen, cd, -bez, t32a);
        t32blen = scaleExpansionZeroElim(cdlen, cd, -beztail, t32b);
        t64alen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64a);
        t32alen = scaleExpansionZeroElim(bdlen, bd, cez, t32a);
        t32blen = scaleExpansionZeroElim(bdlen, bd, ceztail, t32b);
        t64blen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64b);
        t32alen = scaleExpansionZeroElim(bclen, bc, -dez, t32a);
        t32blen = scaleExpansionZeroElim(bclen, bc, -deztail, t32b);
        t64clen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64c);
        t128len = expansionSumZeroElimFast(t64alen, t64a, t64blen, t64b, t128);
        t192len = expansionSumZeroElimFast(t64clen, t64c, t128len, t128, t192);
        
        double[] detx = new double[384];
        double[] detxx = new double[768];
        double[] detxt = new double[384];
        double[] detxxt = new double[768];
        double[] detxtxt = new double[768];
        double[] x1 = new double[1536];
        double[] x2 = new double[2304];
        int xlen = scaleExpansionZeroElim(t192len, t192, aex, detx);
        int xxlen = scaleExpansionZeroElim(xlen, detx, aex, detxx);
        int xtlen = scaleExpansionZeroElim(t192len, t192, aextail, detxt);
        int xxtlen = scaleExpansionZeroElim(xtlen, detxt, aex, detxxt);
        for (int i = 0; i < xxtlen; ++i) {
            detxxt[i] *= 2.0;
        }
        int xtxtlen = scaleExpansionZeroElim(xtlen, detxt, aextail, detxtxt);
        int x1len = expansionSumZeroElimFast(xxlen, detxx, xxtlen, detxxt, x1);
        int x2len = expansionSumZeroElimFast(x1len, x1, xtxtlen, detxtxt, x2);
        
        double[] dety = new double[384];
        double[] detyy = new double[768];
        double[] detyt = new double[384];
        double[] detyyt = new double[768];
        double[] detytyt = new double[768];
        double[] y1 = new double[1536];
        double[] y2 = new double[2304];
        int ylen = scaleExpansionZeroElim(t192len, t192, aey, dety);
        int yylen = scaleExpansionZeroElim(ylen, dety, aey, detyy);
        int ytlen = scaleExpansionZeroElim(t192len, t192, aeytail, detyt);
        int yytlen = scaleExpansionZeroElim(ytlen, detyt, aey, detyyt);
        for (int i = 0; i < yytlen; ++i) {
            detyyt[i] *= 2.0;
        }
        int ytytlen = scaleExpansionZeroElim(ytlen, detyt, aeytail, detytyt);
        int y1len = expansionSumZeroElimFast(yylen, detyy, yytlen, detyyt, y1);
        int y2len = expansionSumZeroElimFast(y1len, y1, ytytlen, detytyt, y2);
        
        double[] detz = new double[384];
        double[] detzz = new double[768];
        double[] detzt = new double[384];
        double[] detzzt = new double[768];
        double[] detztzt = new double[768];
        double[] z1 = new double[1536];
        double[] z2 = new double[2304];
        int zlen = scaleExpansionZeroElim(t192len, t192, aez, detz);
        int zzlen = scaleExpansionZeroElim(zlen, detz, aez, detzz);
        int ztlen = scaleExpansionZeroElim(t192len, t192, aeztail, detzt);
        int zztlen = scaleExpansionZeroElim(ztlen, detzt, aez, detzzt);
        for (int i = 0; i < zztlen; ++i) {
            detzzt[i] *= 2.0;
        }
        int ztztlen = scaleExpansionZeroElim(ztlen, detzt, aeztail, detztzt);
        int z1len = expansionSumZeroElimFast(zzlen, detzz, zztlen, detzzt, z1);
        int z2len = expansionSumZeroElimFast(z1len, z1, ztztlen, detztzt, z2);
        
        double[] detxy = new double[4608];
        double[] adet = new double[6912];
        double[] bdet = new double[6912];
        double[] cdet = new double[6912];
        double[] ddet = new double[6912];
        int xylen = expansionSumZeroElimFast(x2len, x2, y2len, y2, detxy);
        int alen = expansionSumZeroElimFast(z2len, z2, xylen, detxy, adet);
        
        t32alen = scaleExpansionZeroElim(dalen, da, cez, t32a);
        t32blen = scaleExpansionZeroElim(dalen, da, ceztail, t32b);
        t64alen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64a);
        t32alen = scaleExpansionZeroElim(aclen, ac, dez, t32a);
        t32blen = scaleExpansionZeroElim(aclen, ac, deztail, t32b);
        t64blen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64b);
        t32alen = scaleExpansionZeroElim(cdlen, cd, aez, t32a);
        t32blen = scaleExpansionZeroElim(cdlen, cd, aeztail, t32b);
        t64clen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64c);
        t128len = expansionSumZeroElimFast(t64alen, t64a, t64blen, t64b, t128);
        t192len = expansionSumZeroElimFast(t64clen, t64c, t128len, t128, t192);
        xlen = scaleExpansionZeroElim(t192len, t192, bex, detx);
        xxlen = scaleExpansionZeroElim(xlen, detx, bex, detxx);
        xtlen = scaleExpansionZeroElim(t192len, t192, bextail, detxt);
        xxtlen = scaleExpansionZeroElim(xtlen, detxt, bex, detxxt);
        for (int i = 0; i < xxtlen; ++i) {
            detxxt[i] *= 2.0;
        }
        xtxtlen = scaleExpansionZeroElim(xtlen, detxt, bextail, detxtxt);
        x1len = expansionSumZeroElimFast(xxlen, detxx, xxtlen, detxxt, x1);
        x2len = expansionSumZeroElimFast(x1len, x1, xtxtlen, detxtxt, x2);
        ylen = scaleExpansionZeroElim(t192len, t192, bey, dety);
        yylen = scaleExpansionZeroElim(ylen, dety, bey, detyy);
        ytlen = scaleExpansionZeroElim(t192len, t192, beytail, detyt);
        yytlen = scaleExpansionZeroElim(ytlen, detyt, bey, detyyt);
        for (int i = 0; i < yytlen; ++i) {
            detyyt[i] *= 2.0;
        }
        ytytlen = scaleExpansionZeroElim(ytlen, detyt, beytail, detytyt);
        y1len = expansionSumZeroElimFast(yylen, detyy, yytlen, detyyt, y1);
        y2len = expansionSumZeroElimFast(y1len, y1, ytytlen, detytyt, y2);
        zlen = scaleExpansionZeroElim(t192len, t192, bez, detz);
        zzlen = scaleExpansionZeroElim(zlen, detz, bez, detzz);
        ztlen = scaleExpansionZeroElim(t192len, t192, beztail, detzt);
        zztlen = scaleExpansionZeroElim(ztlen, detzt, bez, detzzt);
        for (int i = 0; i < zztlen; ++i) {
            detzzt[i] *= 2.0;
        }
        ztztlen = scaleExpansionZeroElim(ztlen, detzt, beztail, detztzt);
        z1len = expansionSumZeroElimFast(zzlen, detzz, zztlen, detzzt, z1);
        z2len = expansionSumZeroElimFast(z1len, z1, ztztlen, detztzt, z2);
        xylen = expansionSumZeroElimFast(x2len, x2, y2len, y2, detxy);
        int blen = expansionSumZeroElimFast(z2len, z2, xylen, detxy, bdet);
        
        t32alen = scaleExpansionZeroElim(ablen, ab, -dez, t32a);
        t32blen = scaleExpansionZeroElim(ablen, ab, -deztail, t32b);
        t64alen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64a);
        t32alen = scaleExpansionZeroElim(bdlen, bd, -aez, t32a);
        t32blen = scaleExpansionZeroElim(bdlen, bd, -aeztail, t32b);
        t64blen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64b);
        t32alen = scaleExpansionZeroElim(dalen, da, -bez, t32a);
        t32blen = scaleExpansionZeroElim(dalen, da, -beztail, t32b);
        t64clen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64c);
        t128len = expansionSumZeroElimFast(t64alen, t64a, t64blen, t64b, t128);
        t192len = expansionSumZeroElimFast(t64clen, t64c, t128len, t128, t192);
        xlen = scaleExpansionZeroElim(t192len, t192, cex, detx);
        xxlen = scaleExpansionZeroElim(xlen, detx, cex, detxx);
        xtlen = scaleExpansionZeroElim(t192len, t192, cextail, detxt);
        xxtlen = scaleExpansionZeroElim(xtlen, detxt, cex, detxxt);
        for (int i = 0; i < xxtlen; ++i) {
            detxxt[i] *= 2.0;
        }
        xtxtlen = scaleExpansionZeroElim(xtlen, detxt, cextail, detxtxt);
        x1len = expansionSumZeroElimFast(xxlen, detxx, xxtlen, detxxt, x1);
        x2len = expansionSumZeroElimFast(x1len, x1, xtxtlen, detxtxt, x2);
        ylen = scaleExpansionZeroElim(t192len, t192, cey, dety);
        yylen = scaleExpansionZeroElim(ylen, dety, cey, detyy);
        ytlen = scaleExpansionZeroElim(t192len, t192, ceytail, detyt);
        yytlen = scaleExpansionZeroElim(ytlen, detyt, cey, detyyt);
        for (int i = 0; i < yytlen; ++i) {
            detyyt[i] *= 2.0;
        }
        ytytlen = scaleExpansionZeroElim(ytlen, detyt, ceytail, detytyt);
        y1len = expansionSumZeroElimFast(yylen, detyy, yytlen, detyyt, y1);
        y2len = expansionSumZeroElimFast(y1len, y1, ytytlen, detytyt, y2);
        zlen = scaleExpansionZeroElim(t192len, t192, cez, detz);
        zzlen = scaleExpansionZeroElim(zlen, detz, cez, detzz);
        ztlen = scaleExpansionZeroElim(t192len, t192, ceztail, detzt);
        zztlen = scaleExpansionZeroElim(ztlen, detzt, cez, detzzt);
        for (int i = 0; i < zztlen; ++i) {
            detzzt[i] *= 2.0;
        }
        ztztlen = scaleExpansionZeroElim(ztlen, detzt, ceztail, detztzt);
        z1len = expansionSumZeroElimFast(zzlen, detzz, zztlen, detzzt, z1);
        z2len = expansionSumZeroElimFast(z1len, z1, ztztlen, detztzt, z2);
        xylen = expansionSumZeroElimFast(x2len, x2, y2len, y2, detxy);
        int clen = expansionSumZeroElimFast(z2len, z2, xylen, detxy, cdet);
        
        t32alen = scaleExpansionZeroElim(bclen, bc, aez, t32a);
        t32blen = scaleExpansionZeroElim(bclen, bc, aeztail, t32b);
        t64alen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64a);
        t32alen = scaleExpansionZeroElim(aclen, ac, -bez, t32a);
        t32blen = scaleExpansionZeroElim(aclen, ac, -beztail, t32b);
        t64blen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64b);
        t32alen = scaleExpansionZeroElim(ablen, ab, cez, t32a);
        t32blen = scaleExpansionZeroElim(ablen, ab, ceztail, t32b);
        t64clen = expansionSumZeroElimFast(t32alen, t32a, t32blen, t32b, t64c);
        t128len = expansionSumZeroElimFast(t64alen, t64a, t64blen, t64b, t128);
        t192len = expansionSumZeroElimFast(t64clen, t64c, t128len, t128, t192);
        xlen = scaleExpansionZeroElim(t192len, t192, dex, detx);
        xxlen = scaleExpansionZeroElim(xlen, detx, dex, detxx);
        xtlen = scaleExpansionZeroElim(t192len, t192, dextail, detxt);
        xxtlen = scaleExpansionZeroElim(xtlen, detxt, dex, detxxt);
        for (int i = 0; i < xxtlen; ++i) {
            detxxt[i] *= 2.0;
        }
        xtxtlen = scaleExpansionZeroElim(xtlen, detxt, dextail, detxtxt);
        x1len = expansionSumZeroElimFast(xxlen, detxx, xxtlen, detxxt, x1);
        x2len = expansionSumZeroElimFast(x1len, x1, xtxtlen, detxtxt, x2);
        ylen = scaleExpansionZeroElim(t192len, t192, dey, dety);
        yylen = scaleExpansionZeroElim(ylen, dety, dey, detyy);
        ytlen = scaleExpansionZeroElim(t192len, t192, deytail, detyt);
        yytlen = scaleExpansionZeroElim(ytlen, detyt, dey, detyyt);
        for (int i = 0; i < yytlen; ++i) {
            detyyt[i] *= 2.0;
        }
        ytytlen = scaleExpansionZeroElim(ytlen, detyt, deytail, detytyt);
        y1len = expansionSumZeroElimFast(yylen, detyy, yytlen, detyyt, y1);
        y2len = expansionSumZeroElimFast(y1len, y1, ytytlen, detytyt, y2);
        zlen = scaleExpansionZeroElim(t192len, t192, dez, detz);
        zzlen = scaleExpansionZeroElim(zlen, detz, dez, detzz);
        ztlen = scaleExpansionZeroElim(t192len, t192, deztail, detzt);
        zztlen = scaleExpansionZeroElim(ztlen, detzt, dez, detzzt);
        for (int i = 0; i < zztlen; ++i) {
            detzzt[i] *= 2.0;
        }
        ztztlen = scaleExpansionZeroElim(ztlen, detzt, deztail, detztzt);
        z1len = expansionSumZeroElimFast(zzlen, detzz, zztlen, detzzt, z1);
        z2len = expansionSumZeroElimFast(z1len, z1, ztztlen, detztzt, z2);
        xylen = expansionSumZeroElimFast(x2len, x2, y2len, y2, detxy);
        int dlen = expansionSumZeroElimFast(z2len, z2, xylen, detxy, ddet);
        
        double[] abdet = new double[13824];
        double[] cddet = new double[13824];
        double[] det = new double[27648];
        ablen = expansionSumZeroElimFast(alen, adet, blen, bdet, abdet);
        cdlen = expansionSumZeroElimFast(clen, cdet, dlen, ddet, cddet);
        int detlen = expansionSumZeroElimFast(ablen, abdet, cdlen, cddet, det);
        
        return det[detlen - 1];
    }
    
    private static double leftOfPlaneExact(double xa, double ya, double za, double xb, double yb, double zb, double xc,
                                           double yc, double zc, double xd, double yd, double zd) {
        Two t = new Two();
        twoDiff(xa, xd, t);
        double adx = t.x;
        double adxtail = t.y;
        twoDiff(ya, yd, t);
        double ady = t.x;
        double adytail = t.y;
        twoDiff(za, zd, t);
        double adz = t.x;
        double adztail = t.y;
        twoDiff(xb, xd, t);
        double bdx = t.x;
        double bdxtail = t.y;
        twoDiff(yb, yd, t);
        double bdy = t.x;
        double bdytail = t.y;
        twoDiff(zb, zd, t);
        double bdz = t.x;
        double bdztail = t.y;
        twoDiff(xc, xd, t);
        double cdx = t.x;
        double cdxtail = t.y;
        twoDiff(yc, yd, t);
        double cdy = t.x;
        double cdytail = t.y;
        twoDiff(zc, zd, t);
        double cdz = t.x;
        double cdztail = t.y;
        
        double[] axby = new double[8];
        twoTwoProduct(adx, adxtail, bdy, bdytail, axby);
        double negate = -ady;
        double negatetail = -adytail;
        double[] bxay = new double[8];
        twoTwoProduct(bdx, bdxtail, negate, negatetail, bxay);
        
        double[] bxcy = new double[8];
        twoTwoProduct(bdx, bdxtail, cdy, cdytail, bxcy);
        negate = -bdy;
        negatetail = -bdytail;
        double[] cxby = new double[8];
        twoTwoProduct(cdx, cdxtail, negate, negatetail, cxby);
        
        double[] cxay = new double[8];
        twoTwoProduct(cdx, cdxtail, ady, adytail, cxay);
        negate = -cdy;
        negatetail = -cdytail;
        double[] axcy = new double[8];
        twoTwoProduct(adx, adxtail, negate, negatetail, axcy);
        
        double[] t16 = new double[16];
        double[] t32 = new double[32];
        double[] t32t = new double[32];
        int t16len, t32len, t32tlen;
        
        t16len = expansionSumZeroElimFast(8, bxcy, 8, cxby, t16);
        t32len = scaleExpansionZeroElim(t16len, t16, adz, t32);
        t32tlen = scaleExpansionZeroElim(t16len, t16, adztail, t32t);
        double[] adet = new double[64];
        int alen = expansionSumZeroElimFast(t32len, t32, t32tlen, t32t, adet);
        
        t16len = expansionSumZeroElimFast(8, cxay, 8, axcy, t16);
        t32len = scaleExpansionZeroElim(t16len, t16, bdz, t32);
        t32tlen = scaleExpansionZeroElim(t16len, t16, bdztail, t32t);
        double[] bdet = new double[64];
        int blen = expansionSumZeroElimFast(t32len, t32, t32tlen, t32t, bdet);
        
        t16len = expansionSumZeroElimFast(8, axby, 8, bxay, t16);
        t32len = scaleExpansionZeroElim(t16len, t16, cdz, t32);
        t32tlen = scaleExpansionZeroElim(t16len, t16, cdztail, t32t);
        double[] cdet = new double[64];
        int clen = expansionSumZeroElimFast(t32len, t32, t32tlen, t32t, cdet);
        
        double[] abdet = new double[128];
        int ablen = expansionSumZeroElimFast(alen, adet, blen, bdet, abdet);
        double[] det = new double[192];
        int detlen = expansionSumZeroElimFast(ablen, abdet, clen, cdet, det);
        
        return det[detlen - 1];
    }
    
    private static int scaleExpansionZeroElim(int elen, double[] e, double b, double[] h) {
        Two t = new Two();
        split(b, t);
        double bhi = t.x;
        double blo = t.y;
        twoProduct1Presplit(e[0], b, bhi, blo, t);
        double q = t.x;
        double hh = t.y;
        int hindex = 0;
        if (hh != 0) {
            h[hindex++] = hh;
        }
        for (int eindex = 1; eindex < elen; ++eindex) {
            double enow = e[eindex];
            twoProduct1Presplit(enow, b, bhi, blo, t);
            double product1 = t.x;
            double product0 = t.y;
            twoSum(q, product0, t);
            double sum = t.x;
            hh = t.y;
            if (hh != 0) {
                h[hindex++] = hh;
            }
            twoSumFast(product1, sum, t);
            q = t.x;
            hh = t.y;
            if (hh != 0) {
                h[hindex++] = hh;
            }
        }
        if (q != 0.0 || hindex == 0) {
            h[hindex++] = q;
        }
        return hindex;
    }
    
    private static int expansionSumZeroElimFast(int elen, double[] e, int flen, double[] f, double[] h) {
        double q, qnew, hh;
        Two t = new Two();
        double enow = e[0];
        double fnow = f[0];
        int eindex = 0;
        int findex = 0;
        if (fnow > enow == fnow > -enow) {
            q = enow;
            ++eindex;
            if (eindex < elen) {
                enow = e[eindex];
            }
        } else {
            q = fnow;
            ++findex;
            if (findex < flen) {
                fnow = f[findex];
            }
        }
        int hindex = 0;
        if (eindex < elen && findex < flen) {
            if (fnow > enow == fnow > -enow) {
                twoSumFast(enow, q, t);
                qnew = t.x;
                hh = t.y;
                ++eindex;
                if (eindex < elen) {
                    enow = e[eindex];
                }
            } else {
                twoSumFast(fnow, q, t);
                qnew = t.x;
                hh = t.y;
                ++findex;
                if (findex < flen) {
                    fnow = f[findex];
                }
            }
            q = qnew;
            if (hh != 0.0) {
                h[hindex++] = hh;
            }
            while (eindex < elen && findex < flen) {
                if (fnow > enow == fnow > -enow) {
                    twoSum(q, enow, t);
                    qnew = t.x;
                    hh = t.y;
                    ++eindex;
                    if (eindex < elen) {
                        enow = e[eindex];
                    }
                } else {
                    twoSum(q, fnow, t);
                    qnew = t.x;
                    hh = t.y;
                    ++findex;
                    if (findex < flen) {
                        fnow = f[findex];
                    }
                }
                q = qnew;
                if (hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while (eindex < elen) {
            twoSum(q, enow, t);
            qnew = t.x;
            hh = t.y;
            ++eindex;
            if (eindex < elen) {
                enow = e[eindex];
            }
            q = qnew;
            if (hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while (findex < flen) {
            twoSum(q, fnow, t);
            qnew = t.x;
            hh = t.y;
            ++findex;
            if (findex < flen) {
                fnow = f[findex];
            }
            q = qnew;
            if (hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if (q != 0.0 || hindex == 0) {
            h[hindex++] = q;
        }
        return hindex;
    }
    
    private static void twoTwoProduct(double a1, double a0, double b1, double b0, double[] x) {
        double u0, u1, u2, ui, uj, uk, ul, um, un;
        Two t = new Two();
        split(a0, t);
        double a0hi = t.x;
        double a0lo = t.y;
        split(b0, t);
        double b0hi = t.x;
        double b0lo = t.y;
        twoProduct2Presplit(a0, a0hi, a0lo, b0, b0hi, b0lo, t);
        ui = t.x;
        x[0] = t.y;
        split(a1, t);
        double a1hi = t.x;
        double a1lo = t.y;
        twoProduct2Presplit(a1, a1hi, a1lo, b0, b0hi, b0lo, t);
        uj = t.x;
        u0 = t.y;
        twoSum(ui, u0, t);
        uk = t.x;
        u1 = t.y;
        twoSumFast(uj, uk, t);
        ul = t.x;
        u2 = t.y;
        split(b1, t);
        double b1hi = t.x;
        double b1lo = t.y;
        twoProduct2Presplit(a0, a0hi, a0lo, b1, b1hi, b1lo, t);
        ui = t.x;
        u0 = t.y;
        twoSum(u1, u0, t);
        uk = t.x;
        x[1] = t.y;
        twoSum(u2, uk, t);
        uj = t.x;
        u1 = t.y;
        twoSum(ul, uj, t);
        um = t.x;
        u2 = t.y;
        twoProduct2Presplit(a1, a1hi, a1lo, b1, b1hi, b1lo, t);
        uj = t.x;
        u0 = t.y;
        twoSum(ui, u0, t);
        un = t.x;
        u0 = t.y;
        twoSum(u1, u0, t);
        ui = t.x;
        x[2] = t.y;
        twoSum(u2, ui, t);
        uk = t.x;
        u1 = t.y;
        twoSum(um, uk, t);
        ul = t.x;
        u2 = t.y;
        twoSum(uj, un, t);
        uk = t.x;
        u0 = t.y;
        twoSum(u1, u0, t);
        uj = t.x;
        x[3] = t.y;
        twoSum(u2, uj, t);
        ui = t.x;
        u1 = t.y;
        twoSum(ul, ui, t);
        um = t.x;
        u2 = t.y;
        twoSum(u1, uk, t);
        ui = t.x;
        x[4] = t.y;
        twoSum(u2, ui, t);
        uk = t.x;
        x[5] = t.y;
        twoSum(um, uk, t);
        x[7] = t.x;
        x[6] = t.y;
    }
    
    private static void split(double a, Two t) {
        double c = SPLITTER * a;
        double abig = c - a;
        t.x = c - abig;
        t.y = a - t.x;
    }
    private static void twoDiff(double a, double b, Two t) {
        double x = a - b;
        double bvirt = a - x;
        double avirt = x + bvirt;
        double bround = bvirt - b;
        double around = a - avirt;
        t.x = x;
        t.y = around + bround;
    }
    private static void twoSum(double a, double b, Two t) {
        double x = a + b;
        double bvirt = x - a;
        double avirt = x - bvirt;
        double bround = b - bvirt;
        double around = a - avirt;
        t.x = x;
        t.y = around + bround;
    }
    private static void twoSumFast(double a, double b, Two t) {
        double x = a + b;
        double bvirt = x - a;
        t.x = x;
        t.y = b - bvirt;
    }
    private static void twoProduct1Presplit(double a, double b, double bhi, double blo, Two t) {
        split(a, t);
        double ahi = t.x;
        double alo = t.y;
        t.x = a * b;
        double err1 = t.x - ahi * bhi;
        double err2 = err1 - alo * bhi;
        double err3 = err2 - ahi * blo;
        t.y = alo * blo - err3;
    }
    private static void twoProduct2Presplit(double a, double ahi, double alo, double b, double bhi, double blo, Two t) {
        t.x = a * b;
        double err1 = t.x - ahi * bhi;
        double err2 = err1 - alo * bhi;
        double err3 = err2 - ahi * blo;
        t.y = alo * blo - err3;
    }
}
