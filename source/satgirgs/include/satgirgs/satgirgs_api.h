
#ifndef SATGIRGS_API_H
#define SATGIRGS_API_H

#ifdef SATGIRGS_STATIC_DEFINE
#  define SATGIRGS_API
#  define SATGIRGS_NO_EXPORT
#else
#  ifndef SATGIRGS_API
#    ifdef satgirgs_EXPORTS
        /* We are building this library */
#      define SATGIRGS_API __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define SATGIRGS_API __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef SATGIRGS_NO_EXPORT
#    define SATGIRGS_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef SATGIRGS_DEPRECATED
#  define SATGIRGS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef SATGIRGS_DEPRECATED_EXPORT
#  define SATGIRGS_DEPRECATED_EXPORT SATGIRGS_API SATGIRGS_DEPRECATED
#endif

#ifndef SATGIRGS_DEPRECATED_NO_EXPORT
#  define SATGIRGS_DEPRECATED_NO_EXPORT SATGIRGS_NO_EXPORT SATGIRGS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SATGIRGS_NO_DEPRECATED
#    define SATGIRGS_NO_DEPRECATED
#  endif
#endif

#endif /* SATGIRGS_API_H */
