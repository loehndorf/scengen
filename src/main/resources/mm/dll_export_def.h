// ----------------------------------------------------------------
// macros for cross-platform dll-exporting
//
// source: http://gcc.gnu.org/wiki/Visibility

#ifndef DLL_EXPORT_DEF_H
#define DLL_EXPORT_DEF_H

#if defined NO_DLL_DEFS
	// use this if we are just building an executable that does not use DLLs
	#define DLL_PUBLIC
	#define DLL_LOCAL
#else
	#if defined _WIN32
		#ifdef BUILDING_DLL
			#ifdef __GNUC__
				#define DLL_PUBLIC __attribute__((dllexport))
			#else
				// Note: actually gcc seems to also supports this syntax.
				#define DLL_PUBLIC __declspec(dllexport)
			#endif
		#else
			#ifdef __GNUC__
				#define DLL_PUBLIC __attribute__((dllimport))
			#else
				// Note: actually gcc seems to also supports this syntax.
				#define DLL_PUBLIC __declspec(dllimport)
			#endif
		#endif
		#define DLL_LOCAL
	#else
		#if __GNUC__ >= 4
			#define DLL_PUBLIC __attribute__ ((visibility("default")))
			#define DLL_LOCAL  __attribute__ ((visibility("hidden")))
		#else
			#define DLL_PUBLIC
			#define DLL_LOCAL
		#endif
	#endif
#endif

#endif
