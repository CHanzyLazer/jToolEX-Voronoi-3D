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
package jtool.code.collection;

import jtool.code.functional.IOperator1;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Unmodifiable;

import java.util.*;
import java.util.stream.Stream;

/**
 * 获取抽象容器的类，这里获取的结果统一为引用
 * @author CHanzy
 */
public class AbstractCollections {
    private AbstractCollections() {}
    
    /**
     * map {@code Iterable<T> to Iterable<R>} like {@link Stream}.map
     * @author CHanzy
     */
    public static <R, T> @Unmodifiable Iterable<R> map(final Iterable<T> aIterable, final IOperator1<? extends R, ? super T> aOpt) {
        return () -> map(aIterable.iterator(), aOpt);
    }
    public static <R, T> @Unmodifiable Iterator<R> map(final Iterator<T> aIterator, final IOperator1<? extends R, ? super T> aOpt) {
        return new Iterator<R>() {
            final Iterator<T> mIt = aIterator;
            @Override public boolean hasNext() {
                return mIt.hasNext();
            }
            @Override public R next() {
                if (hasNext()) {
                    return aOpt.cal(mIt.next());
                } else {
                    throw new NoSuchElementException();
                }
            }
        };
    }
    public static <R, T> @Unmodifiable Collection<R> map(final Collection<T> aCollection, final IOperator1<? extends R, ? super T> aOpt) {
        return new AbstractCollection<R>() {
            @Override public @NotNull Iterator<R> iterator() {return map(aCollection.iterator(), aOpt);}
            @Override public int size() {return aCollection.size();}
        };
    }
    public static <R, T> @Unmodifiable List<R> map(final List<T> aList, final IOperator1<? extends R, ? super T> aOpt) {
        return new AbstractRandomAccessList<R>() {
            @Override public R get(int index) {return aOpt.cal(aList.get(index));}
            @Override public int size() {return aList.size();}
            @Override public @NotNull Iterator<R> iterator() {return map(aList.iterator(), aOpt);}
        };
    }
    public static <R, T> @Unmodifiable List<R> map(final T[] aArray, final IOperator1<? extends R, ? super T> aOpt) {
        return new AbstractRandomAccessList<R>() {
            @Override public R get(int index) {return aOpt.cal(aArray[index]);}
            @Override public int size() {return aArray.length;}
        };
    }
}
