#ifdef DEFINE_GLOBAL
#     define EXTERN
#else
#     define EXTERN extern
#endif

EXTERN int        ndim, nvol, *lsize, **nn;

#undef EXTERN
