#include "core/atomic.h"

#include <atomic>

struct atomic_int {
	std::atomic<int> count {1};
};

extern "C" {
atomic_int* atomic_int_alloc() {
	return new atomic_int;
}

int atomic_int_decref(atomic_int* a) {
	int n = --a->count;
	if (n == 0)
		delete a;
	return n;
}

int atomic_int_incref(atomic_int* a) {
	return ++a->count;
}
};
