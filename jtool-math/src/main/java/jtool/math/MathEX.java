/**
 * Copyright (C) 2023 CHanzy/CHanzyLazer. All rights reserved.
 *
 * This file is part of jtool
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package jtool.math;

import jtool.atom.IXYZ;
import jtool.atom.XYZ;
import jtoolex.algorithm.Geometry;

/**
 * @author CHanzy
 * <p> Extended mathematical methods </p>
 * <p> The method of using internal Thread Pool is not thread safe when {@code nThreads > 1} </p>
 */
public class MathEX {
    /** 使用外套一层的写法可以在外部使用是替换底层的几何算法，使用优化的版本 */
    public static class Graph {
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
            return Geometry.leftOfPlane(
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
            return Geometry.inSphere(
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
            Geometry.centerSphere(
                  aA.x(), aA.y(), aA.z()
                , aB.x(), aB.y(), aB.z()
                , aC.x(), aC.y(), aC.z()
                , aD.x(), aD.y(), aD.z()
                , rCenter);
            return rCenter;
        }
    }
}
