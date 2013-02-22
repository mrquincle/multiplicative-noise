/*
 * helix.h
 *
 * New type of data structure. I've never seen it before myself. And I'm pretty
 * proud on it. :-) It can integrate 1-dimensional arrays very fast!
 *
 * This software is published under the GNU Lesser General Public license (LGPL).
 *
 * It is not possible to add usage restrictions to an open-source license. Nevertheless,
 * we personally strongly object against this software used by the military, in the
 * bio-industry, for animal experimentation, or anything that violates the Universal
 * Declaration of Human Rights.
 *
 * Copyright © 2011 Anne van Rossum <anne@almende.com>
 *
 * @author 	Anne C. van Rossum
 * @date	Sep 20, 2011
 * @project	Replicator FP7
 * @company	Almende B.V.
 * @license LGPL v3 or newer
 * @case	Self-organised criticality
 */


#ifndef HELIX_H_
#define HELIX_H_

#include <stddef.h>

/**
 * Specific data structure very convenient for 1-dimensional integration:
 * 1.) No copying around from left <- center <- right but using two arrays
 * 2.) Very fast circular implementation by updating boundary values only in update()
 * 3.) Very fast left() and right() implementation by dedicated shifted array indices.

 * 1.) There are two arrays instead of one (memory footprint is *2). Integration is
 * over old values. Instead of some "sliding window" we just calculate the new value
 * given some neighbours (left and right) and put it in another array. When all values
 * are calculated, the user is responsible for calling "update()" which internally
 * swaps the source and target arrays.
 * 2.) To make it a circular array, no linked list construction is created (with get()
 * overhead), but an array that is two elements larger than the original. The first and
 * last values are copied on update(), so left and right make also sense for first and
 * last item.
 * 3.) There is a "me_right" array that is shifted to the right, so index "i" points
 * directly to the item at the right (without incrementing anything). Likewise with
 * "me_raw" for the item at the left.
 */
template <typename T>
class helix {
public:
	helix(int len) {
		N = len+2;
		// the data vectors are circular (last item is also placed at beginning)
		u = new T[N];
		v = new T[N];
		// we toggle between u and v
		me_raw = u;
		you_raw = v;
		// we need shifted vectors by one for get() and set()
		me = &me_raw[1];
		you = &you_raw[1];

		me_right = &me_raw[2];
	}

	~helix() {
		me = you = me_raw = you_raw = NULL;
		delete [] u;
		delete [] v;
	}

	void update() {
		// update boundaries
		you_raw[0] = you_raw[N-2];
		you_raw[N-1] = you_raw[1];
		// toggle me and you
		T* temp = me_raw;
		me_raw = you_raw;
		you_raw = temp;
		// give pointers to shifted vectors
		me = &me_raw[1];
		you = &you_raw[1];

		me_right = &me_raw[2];
	}

	inline T get(int i) {
		return me[i];
	}

//	inline T raw(int i) {
//		return me_raw[i];
//	}

	// preferably I would like to overload operator with and without lvalue
	// but that is not possible, so just use get() and set()

//	T operator[](const int i) { return me[i]; }
//	T& operator[](const int i) { return you[i]; }

	inline void set(T value, int i) {
		you[i] = value;
	}

	//! Get left item
	inline T left(int i) {
		return me_raw[i];
	}

	//! Get right item
	inline T right(int i) {
		return me_right[i];
	}

private:
	T *u;

	T *v;

	int N;

	T *me;

	T *you;

	T *me_raw; // is also me_left

	T *you_raw;

	T *me_right;
};

#endif /* HELIX_H_ */
