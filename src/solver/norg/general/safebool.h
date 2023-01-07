/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _SAFEBOOL_H_
#define _SAFEBOOL_H_

/* The Safe Bool Idiom */

/*
There are two plausible solutions: Using a base class with a virtual function for the actual logic, 
or a base class that knows which function to call on the derived class. As virtual functions come 
at a cost (especially if the class you're augmenting with Boolean tests doesn't contain any other 
virtual functions), I add support for both versions below:
*/

class safe_bool_base {
protected:
	typedef void (safe_bool_base::*bool_type)() const;
	void this_type_does_not_support_comparisons() const {}
	safe_bool_base() {}
	safe_bool_base(const safe_bool_base &) {}
	safe_bool_base &operator=(const safe_bool_base &) {return *this;}
	~safe_bool_base() {}
  };

template <typename T>
class safe_bool: public safe_bool_base {
public:
	operator bool_type() const {
		return (static_cast<const T *>(this))->boolean_test()
			? &safe_bool::this_type_does_not_support_comparisons : 0;
	}
protected:
	~safe_bool() {}
};

template <>
class safe_bool<void>: public safe_bool_base {
public:
	operator bool_type() const {
		return boolean_test() == true ?
			&safe_bool::this_type_does_not_support_comparisons : 0;
	}
protected:
	virtual bool boolean_test() const = 0;
	virtual ~safe_bool() {}
};

template <typename T, typename U> 
bool operator==(const safe_bool<T> &lhs, const safe_bool<U> &rhs) {
	lhs.this_type_does_not_support_comparisons();
	return false;
}

template <typename T,typename U> 
bool operator!=(const safe_bool<T> &lhs, const safe_bool<U> &rhs) {
	lhs.this_type_does_not_support_comparisons();
	return false;	
}

/*
// example 1
class Testable_with_virtual: public safe_bool<void> {
protected:
	bool boolean_test() const {
		// Perform Boolean logic here
	}
};

// example 2 (this one is suggested)
class Testable_without_virtual: public safe_bool<Testable_without_virtual> {
public:
	bool boolean_test() const {
		// Perform Boolean logic here
	}
};

The first class, Testable_with_virtual, derives publicly from safe_bool, 
and implements a virtual function boolean_test¡ªthis function is called 
whenever an instance is tested (as in if (obj) {}, or if (!obj) {}). 
The second class, Testable_without_virtual, also derives publicly from safe_bool, 
and in addition, it passes itself as a template parameter to its base class. 
This little trick--known as the Curiously Recurring Template Pattern--enables 
the base class to downcast (to the derived class) using static_cast and call 
boolean_test with no extra runtime overhead and no virtual function calls. 
Some people may feel that this is a slight misuse of inheritance; while it 
might be argued that an instance of a derived class is-a safe_bool of sorts, 
this is certainly not the intent of this code. However, there is little reason 
to believe that even neophyte programmers will fall into the trap of 
misunderstanding this relationship. The destructors of the safe_bool classes 
are protected to minimize the potential for misuse. But there's still hope for 
the (in my opinion, overly) conscientious object-oriented purist; use private 
inheritance, and make the conversion function public by reintroducing it in 
the correct scope:

// example 3
class Testable_without_virtual: private safe_bool<Testable_without_virtual> {
public:
	using safe_bool<Testable_without_virtual>::operator bool_type;
	bool boolean_test() const {
		return true; // Logic goes here!
	}
};

Matthew Wilson [5] pointed out that the inheritance strategy (using safe_bool as a base class) 
may lead to size penalties on some compilers, specifically, those that do not implement EBO 
(Empty Base Optimization) properly. Although most modern compilers do when it comes to single 
inheritance, there may be a size penalty with multiple inheritance.
*/

#endif /* _SAFEBOOL_H_ */
