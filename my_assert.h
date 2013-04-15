#ifndef MY_ASSERT_H
#define MY_ASSERT_H

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cerrno>
#include<string>
#include<sstream>
#include<iomanip>

inline std::string itos(int i) {
  std::ostringstream strout;
  strout<<i;
  return strout.str();
}

// n : number of significant digits
inline std::string ftos(double f, int n) {
  std::ostringstream strout;
  strout<<std::setprecision(n)<<f;
  return strout.str();
}

inline std::string ctos(char c) {
  return std::string(1, c);
}

inline std::string cstrtos(const char* s) {
  return std::string(s);
}


#define general_assert(expr, errmsg) if (!(expr)) general_report((errmsg), false)
#define general_assert_1(expr, errmsg) if (!(expr)) general_report((errmsg), true)

inline void general_report(const std::string& errmsg, bool putEnter) {
	if (putEnter) printf("\n");
	fprintf(stderr, "%s\n", errmsg.c_str());
	exit(-1);
}

#define pthread_assert(rc, func_name, errmsg) if ((rc) != 0) pthread_report((rc), (func_name), (errmsg))

inline void pthread_report(int rc, const std::string& func_name, const std::string& errmsg) {
	fprintf(stderr, "%s\n", errmsg.c_str());

	if (func_name == "pthread_create") {
		switch(rc) {
		case EAGAIN:
			fprintf(stderr, "Error code: EAGAIN. Insufficient resources to create another thread, or a system-imposed limit on the number of threads was encountered.\n");
			break;
		case EINVAL:
			fprintf(stderr, "Error code: EINVAL. Invalid settings in attr.\n");
			break;
		case EPERM:
			fprintf(stderr, "Error code: EPERM. No permission to set the scheduling policy and parameters specified in attr.\n");
			break;
		default: fprintf(stderr, "Unknown error code: %d.\n", rc);
		}
	} else if (func_name == "pthread_join") {
		switch(rc) {
		case EDEADLK:
			fprintf(stderr, "Error code: EDEADLK. A deadlock was detected (e.g., two threads tried to join with each other); or thread_id specifies the calling thread.\n");
			break;
		case EINVAL:
			fprintf(stderr, "Error code: EINVAL. The implementation has detected that the value specified by thread_id does not refer to a joinable thread.\n");
			break;
		case ESRCH:
			fprintf(stderr, "Error code: ESRCH. No thread with thread_id could be found.\n");
			break;
		default: fprintf(stderr, "Unknown error code: %d.\n", rc);
		}
	} else if (func_name == "pthread_mutex_lock") {
		switch(rc) {
		case EAGAIN:
			fprintf(stderr, "Error code: EAGAIN. The mutex could not be acquired because the maximum number of recursive locks for mutex has been exceeded.\n");
			break;
		case EDEADLK:
			fprintf(stderr, "Error code: EDEADLK. The current thread already owns the mutex.\n");
			break;
		case EINVAL:
			fprintf(stderr, "Error code: EINVAL. The mutex was created with the protocol attribute having the value PTHREAD_PRIO_PROTECT and the calling thread's priority is higher than the mutex's current priority ceiling. Or the value specified by mutex does not refer to an initialized mutex object.\n");
			break;
		default: fprintf(stderr, "Unknown error code: %d.\n", rc);
		}
	} else if (func_name == "pthread_mutex_unlock") {
		switch(rc) {
		case EAGAIN:
			fprintf(stderr, "Error code: EAGAIN. The mutex could not be acquired because the maximum number of recursive locks for mutex has been exceeded.\n");
			break;
		case EINVAL:
			fprintf(stderr, "Error code: EINVAL. The value specified by mutex does not refer to an initialized mutex object.\n");
			break;
		case EPERM:
			fprintf(stderr, "Error code: EPERM. The current thread does not own the mutex.\n");
			break;
		default: fprintf(stderr, "Unknown error code: %d.\n", rc);
		}
	} else {
		fprintf(stderr, "Unknown function name: %s.\n", func_name.c_str());
	}

	exit(-1);
}

#endif
