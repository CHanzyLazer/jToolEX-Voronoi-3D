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
package jtool.atom;


/**
 * {@link IXYZ} 的一般实现，考虑到效率这里可以直接访问内部成员，从而避免多态函数调用的损失
 * @author CHanzy
 */
public final class XYZ extends AbstractXYZ {
    /**
     * Convert IXYZ to XYZ to optimise, result should be read only!
     * @author CHanzy
     */
    public static XYZ toXYZ(IXYZ aXYZ) {
        return (aXYZ instanceof XYZ) ? (XYZ)aXYZ : new XYZ(aXYZ);
    }
    
    public double mX, mY, mZ;
    public XYZ(double aX, double aY, double aZ) {
        mX = aX; mY = aY; mZ = aZ;
    }
    public XYZ(IXYZ aXYZ) {
        mX = aXYZ.x(); mY = aXYZ.y(); mZ = aXYZ.z();
    }
    public XYZ(XYZ aXYZ) {
        mX = aXYZ.mX; mY = aXYZ.mY; mZ = aXYZ.mZ;
    }
    
    /** 批量设置的接口，返回自身方便链式调用 */
    public XYZ setXYZ(double aX, double aY, double aZ) {
        mX = aX; mY = aY; mZ = aZ;
        return this;
    }
    public XYZ setXYZ(IXYZ aXYZ) {
        mX = aXYZ.x(); mY = aXYZ.y(); mZ = aXYZ.z();
        return this;
    }
    public XYZ setXYZ(XYZ aXYZ) {
        mX = aXYZ.mX; mY = aXYZ.mY; mZ = aXYZ.mZ;
        return this;
    }
    
    
    @Override public double x() {return mX;}
    @Override public double y() {return mY;}
    @Override public double z() {return mZ;}
    
    @Override public double prod() {return mX * mY * mZ;}
    @Override public double min() {return Math.min(Math.min(mX, mY), mZ);}
    @Override public double max() {return Math.max(Math.max(mX, mY), mZ);}
    
    @Override public XYZ cross(double aX, double aY, double aZ) {return new XYZ(mY*aZ - aY*mZ, mZ*aX - aZ*mX, mX*aY - aX*mY);}
    public XYZ cross(XYZ aRHS) {return cross(aRHS.mX, aRHS.mY, aRHS.mZ);}
    public void cross2this(IXYZ aRHS) {cross2this(aRHS.x(), aRHS.y(), aRHS.z());}
    public void cross2this(XYZ aRHS) {cross2this(aRHS.mX, aRHS.mY, aRHS.mZ);}
    public void cross2this(double aX, double aY, double aZ) {
        double tX = mX;
        double tY = mY;
        double tZ = mZ;
        mX = tY*aZ - aY*tZ;
        mY = tZ*aX - aZ*tX;
        mZ = tX*aY - aX*tY;
    }
    @Override public XYZ negative() {return new XYZ(-mX, -mY, -mZ);}
    public void negative2this() {mX = -mX; mY = -mY; mZ = -mZ;}
    @Override public double norm() {return Math.sqrt(mX*mX + mY*mY + mZ*mZ);}
    
    @Override public XYZ plus(IXYZ aRHS) {return new XYZ(mX+aRHS.x(), mY+aRHS.y(), mZ+aRHS.z());}
    @Override public XYZ plus(double aX, double aY, double aZ) {return new XYZ(mX+aX, mY+aY, mZ+aZ);}
    @Override public XYZ plus(double aRHS) {return new XYZ(mX+aRHS, mY+aRHS, mZ+aRHS);}
    /** 使用重载而不是 instanceof，即只优化可以在编译期间判断的情况 */
    public XYZ plus(XYZ aRHS) {return new XYZ(mX+aRHS.mX, mY+aRHS.mY, mZ+aRHS.mZ);}
    public void plus2this(XYZ aRHS) {mX += aRHS.mX; mY += aRHS.mY; mZ += aRHS.mZ;}
    public void plus2this(double aRHS) {mX += aRHS; mY += aRHS; mZ += aRHS;}
    
    @Override public XYZ minus(IXYZ aRHS) {return new XYZ(mX-aRHS.x(), mY-aRHS.y(), mZ-aRHS.z());}
    @Override public XYZ minus(double aX, double aY, double aZ) {return new XYZ(mX-aX, mY-aY, mZ-aZ);}
    @Override public XYZ minus(double aRHS) {return new XYZ(mX-aRHS, mY-aRHS, mZ-aRHS);}
    /** 使用重载而不是 instanceof，即只优化可以在编译期间判断的情况 */
    public XYZ minus(XYZ aRHS) {return new XYZ(mX-aRHS.mX, mY-aRHS.mY, mZ-aRHS.mZ);}
    public void minus2this(XYZ aRHS) {mX -= aRHS.mX; mY -= aRHS.mY; mZ -= aRHS.mZ;}
    public void minus2this(double aRHS) {mX -= aRHS; mY -= aRHS; mZ -= aRHS;}
    
    @Override public XYZ multiply(IXYZ aRHS) {return new XYZ(mX*aRHS.x(), mY*aRHS.y(), mZ*aRHS.z());}
    @Override public XYZ multiply(double aX, double aY, double aZ) {return new XYZ(mX*aX, mY*aY, mZ*aZ);}
    @Override public XYZ multiply(double aRHS) {return new XYZ(mX*aRHS, mY*aRHS, mZ*aRHS);}
    /** 使用重载而不是 instanceof，即只优化可以在编译期间判断的情况 */
    public XYZ multiply(XYZ aRHS) {return new XYZ(mX*aRHS.mX, mY*aRHS.mY, mZ*aRHS.mZ);}
    public void multiply2this(XYZ aRHS) {mX *= aRHS.mX; mY *= aRHS.mY; mZ *= aRHS.mZ;}
    public void multiply2this(double aRHS) {mX *= aRHS; mY *= aRHS; mZ *= aRHS;}
    
    @Override public XYZ div(IXYZ aRHS) {return new XYZ(mX/aRHS.x(), mY/aRHS.y(), mZ/aRHS.z());}
    @Override public XYZ div(double aX, double aY, double aZ) {return new XYZ(mX/aX, mY/aY, mZ/aZ);}
    @Override public XYZ div(double aRHS) {return new XYZ(mX/aRHS, mY/aRHS, mZ/aRHS);}
    /** 使用重载而不是 instanceof，即只优化可以在编译期间判断的情况 */
    public XYZ div(XYZ aRHS) {return new XYZ(mX/aRHS.mX, mY/aRHS.mY, mZ/aRHS.mZ);}
    public void div2this(XYZ aRHS) {mX /= aRHS.mX; mY /= aRHS.mY; mZ /= aRHS.mZ;}
    public void div2this(double aRHS) {mX /= aRHS; mY /= aRHS; mZ /= aRHS;}
    
    
    @Override public double distance2(IXYZ aRHS) {
        double tX = mX - aRHS.x();
        double tY = mY - aRHS.y();
        double tZ = mZ - aRHS.z();
        return tX*tX + tY*tY + tZ*tZ;
    }
    @Override public double distance2(double aX, double aY, double aZ) {
        aX -= mX;
        aY -= mY;
        aZ -= mZ;
        return aX*aX + aY*aY + aZ*aZ;
    }
    /** 使用重载而不是 instanceof，即只优化可以在编译期间判断的情况 */
    public double distance2(XYZ aRHS) {
        double tX = mX - aRHS.mX;
        double tY = mY - aRHS.mY;
        double tZ = mZ - aRHS.mZ;
        return tX*tX + tY*tY + tZ*tZ;
    }
    public double distance(XYZ aRHS) {return Math.sqrt(distance2(aRHS));}
}
