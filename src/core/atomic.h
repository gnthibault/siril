#ifndef CORE_ATOMIC_H
#define CORE_ATOMIC_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct atomic_int atomic_int;

extern atomic_int* atomic_int_alloc();
extern void atomic_int_free(atomic_int* a);
extern int atomic_int_decref(atomic_int* a);
extern int atomic_int_incref(atomic_int* a);

#ifdef __cplusplus
};
#endif

#endif /* CORE_ATOMIC_H */
