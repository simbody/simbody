Contributing to Simbody
=======================
**(work in progress)**

Simbody is a community resource and we encourage you to contribute in whatever way you can -- for example: new code, bug fixes, test cases, and examples; documentation improvements and typo fixes; bug reports, feature requests, ideas and discussion topics; and user forum questions and answers. We appreciate contributions and our development team is collaborative and constructive -- don't be shy! 

**Important note:** Simbody is an open source project licensed under extremely flexible terms intended to encourage use by *anyone*, for *any purpose*. We *want* people to use our software and we *don't want* to restrict how they use it or control their behavior. When you make a contribution to the Simbody project, **you are agreeing** to do so under those same terms. The details are [below](#contributor-license-agreement); if you aren't comfortable with those terms, we're still friends but you shouldn't contribute. 

Contents:
- [Ways to Contribute](#ways-to-contribute)
- [Submitting Pull Requests](#submitting-pull-requests)
- [Coding Conventions](#coding-conventions)
- [List of Contributors](#list-of-contributors)
- [Contributor License Agreement](#contributor-license-agreement)

Ways to contribute
------------------
- Ask and/or answer questions on the [Simbody user forum](https://simtk.org/forums/viewforum.php?f=47).
- File bug reports, documentation problems, feature requests, and developer discussion topics using the GitHub [issue tracker](https://github.com/simbody/simbody/issues).
- Submit GitHub pull requests providing fixes for code or documentation.

Submitting pull requests (PRs)
------------------------------
Please don't surprise us with big out-of-the-blue PRs. If you want to work on an existing issue, post a comment to that effect in the issue. That way you can engage in discussion with others to coordinate your work with theirs and avoid duplication, to see what people think of the approach you have planned, and so on. If there is not yet a relevant issue in place, a great way to start is by creating one, and then engaging in an issue conversation.

The main (upstream) repository (repo) for Simbody is the `simbody` repo in the `simbody` organization on GitHub; that is, https://github.com/simbody/simbody. The general idea is that you will work in your own copy of that repo on GitHub (called a *fork*) and then submit a PR requesting to merge a branch in your fork into a branch of the upstream repo.

### Mechanics of submitting a PR
This is a very abbreviated description of the process. If you are a git newbie you will need to look at some of the great GitHub tutorials, starting with GitHub Bootcamp [here](https://github.com).

Below we'll assume your GitHub account is `yourid`.

1. Create your own fork `yourid/simbody` of the `simbody/simbody` repo on GitHub. Use the `Fork` button [here](https://github.com/simbody/simbody).
2. Clone `yourid/simbody` repo onto your local machine.
3. Create a branch `something-feature` for your new feature or `fix-something-issue123` for a bug fix (we're not fussy about branch names; it is just temporary).
4. Commit the new code or documentation to `something-feature` branch
5. Push `something-feature` branch up to `yourid/simbody` fork on GitHub
6. Go to the `simbody/simbody` repo, click Pull Requests, and create a new PR. Specify `simbody/simbody master` as the base (destination) branch and `yourid/simbody something-feature` as the head (source) branch. 
7. Provide a description and reference the corresponding issue(s).
8. Engage in discussion with Simbody maintainers who will review your changes and make comments. 
9. At the end of the PR Conversation page is a notification that our continuous integration (automatic build) system is trying to build your PR. Check to see whether the build succeeds on all platforms, and if not click the `Details` button and fix the problem if you can, or else ask for help.
  
In most cases discussions and build problems will require you to make some changes to your submission. That is very easy to do because a PR is a *reference* to your branch, not a copy. So you just make the changes to the `something-feature` (or whatever) branch on your local clone, and then push those changes back to the same branch on your `yourid/simbody` fork on GitHub. The changes will immediately start building and you can return to discussing them in the same PR.

Eventually your PR will be merged (good) or closed unmerged by a Simbody maintainer, but always after an open discussion.   

Coding Conventions
------------------
**(still working on this section)**

Existing Simbody code does not perfectly follow these guidelines and we appreciate issues pointing out problems, and especially PRs that correct our earlier slip-ups.

### General code layout

#### Keep line width to 80 characters
Line widths should be no longer than **80** characters. The reason for this is that it permits multiple files to be laid out side-by-side during editing, which is *really* useful. At 80 characters you can get three windows on a modest-sized monitor, using a font that is readable even by adults.

It is best to use a “guide line” (a vertical line that marks column 80 while you edit) so that you can see where the limit is. If you are using Visual Studio, there is a very nice Editor Guidelines Extension available [here](https://visualstudiogallery.msdn.microsoft.com/da227a0b-0e31-4a11-8f6b-3a149cf2e459). If you don't have built-in guide lines available, note that the last line of the copyright block at the top of every Simbody source file is 80 characters wide.

Please don't interpret this to mean we like short lines. On the contrary it is nice to see as much code as reasonably possible on the screen at once, so don't stretch out your code vertically unnecessarily; and don't waste horizontal space where it doesn't help readability. Long comment blocks in particular should have lines as wide as possible. Just don't go over 80.

#### Indent by 4 spaces; no tabs in files
There *must not* be any tabs in your code; our build system will reject them. Please be sure that your code editor is set to replace tabs with four spaces. You should never allow tab characters to get into your code. They will look great to you, but in other viewers people will see your code as randomly formatted.

If you use Visual Studio, go to `Tools:Options:Text Editor:C/C++:Tabs`, set `tab size=indent size=4`, and check the `Insert spaces` button. In `vi` or `vim` use `set tabstop=4` and `set expandtab`. Almost any editor has a similar option, and most can help you clean up a file that already has tabs in it. 

### Naming conventions
We do not believe it is helpful to attempt to encode type information into symbol names (for example, beginning pointer names with a `p`). Much of the need for such conventions has passed with the wide availability of IDEs offering language-sensitive code browsing and debugging, such as that provided by Visual Studio or Eclipse. Thus we do not use name prefix characters to provide information that can easily be obtained while browsing code or debugging. We trust programmers to add appropriate conventions in their own code when those conventions are necessary for clarity or convenience, and to explain them in nearby comments.

We prefer consistency with existing precedent to our own conventions whenever appropriate. For example, C++ containers like `std::vector` define standard names like `const_iterator` so if you are building a container intended to be compatible with one of those you should follow the existing precedent rather than use the Simbody convention which would have been `ConstIterator`.

#### Types: classes and structs, typedefs, enumeration types
Names should be nouns with initial cap and cap for each word (camelcase). There should be no underscores.
```cpp
System
StateClient
Vector
Real
```

#### Constants
We reserve the ugly `ALL_CAPS_CONVENTION` for preprocessor macros both to discourage their use and to signal the type-unsafe loophole being employed. In particular, we discourage the use of preprocessor macros for constants and use a different, less violent convention for type-safe constants.

For constants defined within the language, using `enum` or `const`, use an initial cap, and a cap for each subsequent word (same convention as for classes).
```cpp
enum Color {
    Red,
    Blue,
    LightPink
};
static const Color AttentionColor = Red;
```


#### Functions and methods
Names should begin with a verb, start with lower case, then initial cap for each subsequent word.

```cpp
getBodyVelocity()
setDefaultLength()
```

We have some conventional starting verbs and you should use the same ones when they apply:

   verb   | meaning
----------|---------
`get`     | Return a const reference to something that has already been computed.
`set`     | Change the value of some internal quantity; may cause invalidation of dependent computations. 
`upd`     | (update) Return a writable reference to some internal variable. May cause invalidation of dependent computations.
`find`    | Perform a small calculation (e.g., find the distance between two points) and return the result without changing anything else.
`calc`    | (calculate) Perform an expensive calculation and return the result. Does not cause any other changes.
`realize` | Initiate state-dependent computations and cache results internally; no result returned.
`adopt`   | Take over ownership of a passed-in heap-allocated object.

#### Variables (including local, global, static, members)
Use generally descriptive noun phrases, start with lower case, caps for each word, no underscores. Spell things out unless there is a compelling reason to abbreviate. Follow conventional exceptions for mathematics and indexing.
```cpp
fileName
nextAvailableSlot
i,j,k
```
We do not require that you give data members a distinguishing prefix. However, it is often helpful to do so in complicated classes, in which case the prefix should be `m_`, added to names that otherwise follow the above convention. Do not use an initial underscore alone.
```cpp
m_fileName
m_nextAvailableSlot
```
Please do not use any other prefix conventions; many exist and none are widely agreed upon so they are not helpful to a mixed audience of readers.

#### Preprocessor macro (includes both constants and macros with arguments)
All caps, words separated by underscores. When these appear in interfaces they must be prefixed with a distinctive prefix to keep them from colliding with other symbols. Use an all-caps version of the associated name space when possible. The names of all macros from Simbody software are prefixed with `SimTK_`.

```cpp
SimTK_DEBUG
MYMODULE_UGLY_MACRO
```

#### Namespaces
Short, cryptic, low probability of having the same name as someone else’s namespace. We reserve namespaces containing `SimTK` and `Simbody` (in any combination of upper/lowercase) for user-visible Simbody code.
```cpp
std::
SimTK::
```

### Miscellaneous C++ issues

#### Use anonymous namespaces
If you define classes or external functions in C++ source, even if they appear nowhere else, those names will be exported at link time and may conflict with other names. If that's intentional, make sure the names are in the `SimTK` namespace or begin with `SimTK_`. If not, you should surround the declaration with an anonymous namespace:
```cpp
namespace {
    // declarations that are private to this source file
};
``` 
That prevents the symbols from being exported and you can use any names for them that you want.

For functions you can achieve the same thing by declaring them `static` (which you must do if your code is in C) but anonymous namespaces in C++ are much more powerful.

#### Throw and return are not functions

In C++ `throw` and `return` are not functions. It is misleading to enclose their arguments in parentheses. That is, you should write `return x;` not `return(x);`. A parenthesized expression is not treated the same as a function argument list. For example `f(a,b)` and `return(a,b)` mean very different things -- the former is a 2-argument function call; the latter is an invocation of the rarely-used “comma operator”.

#### Always use pre-increment and pre-decrement operators when you have a choice

Operators for both pre-increment (`++i`) and post-increment (`i++`) are available in C++. If you don’t look at the result, they are logically equivalent. For simple types they are physically equivalent too. But for complicated types (like iterators), the pre-increment operators are much cheaper computationally, because they don’t require separate storage for saving the previous result. Therefore you should get in the habit of using pre-increment (or pre-decrement) in all your loops:

```cpp
for (int i; i < limit; ++i) /* <-- YES*/ 
for (int i; i < limit; i++) /* <-- NO */ 
```

This will prevent you from using the wrong operator in the expensive cases, which are not always obvious. Of course in cases where you actually need the pre- or post-value for something, you should use the appropriate operator. 

#### Place pointer and reference symbols with the type

References and pointers create new types. That is `T`, `T*`, and `T&` are three distinct types. You can tell because you can make `typedef`s like this:

```cpp
typedef T  SameAsT; 
typedef T* PointerToT;
typedef T& ReferenceToT;
 
// and then declare
 
SameAsT      t1,      t2;      // both are type T
PointerToT   tptr1,   tptr2;   // both are type T* 
ReferenceToT tref1=a, tref2=b; // both are type T&
```

Therefore you should place the `*` and `&` next to the type, not the variable, because logically they are part of the type. Unfortunately, the C language had a bug in its syntax which has been inherited by C++. A line like `char* a,b;` is treated like `char* a; char b;` rather than `char* a; char* b;`, but if you write `typedef char* CharPtr;` then `CharPtr a,b;` declares both to be pointers. There is no perfect solution because the language is broken. However, there is no problem in argument lists (since each variable has to have its own type). So we recommend that you simply avoid the misleading multiple-declaration form when using pointers or references. Just use separate declarations or a `typedef`. Then always put the `*` and `&` with the type where they belong. Here are right and wrong examples for argument lists:

```cpp
f(int I, string& name, char* something) /* <-- YES*/ 
f(int I, string &name, char *something) /* <-- NO */ 
```

#### Assignment Operators in C++

You should let the compiler automatically generate the copy constructor and copy assignment operator for your classes whenever possible. But sometimes you have to write one. Here is the basic template for copy assignment:

```cpp
MyClass& operator=(const MyClass& source) {
    if (&source != this) {
       // copy stuff from source to this
    }
    return *this;
} 
```

A common mistake is to leave out the `if`. Since the “copy stuff” part often begins by deleting the contents of “this”, a self assignment like a=a will fail without those lines; that is always supposed to work (and does for all the built-in and standard library types). Of course no one intentionally does that kind of assignment, but they occur anyway in general code since you don’t always know where the source comes from.

If the “copy stuff” part consists only of assignments that work for self assignment, then you can get away without the test, but unless you’ve thought it through carefully you should just get in the habit of putting in the test.


### Use English and a spell checker
While we do encourage international contributors and users, there is much to be gained by choosing a single language, and English is the obvious choice for us. Any submitted code must be understandable by English speakers, with English comments, error messages, and documentation. As one practical consequence, this means Simbody code can use `char` rather than `wchar_t` (wide characters) and embedded English text to be displayed at run time is acceptable and need not be sequestered into separate files to facilitate translation. Simbody contributors thus do not need to be familiar with techniques for international programming.

Please use correct spelling and punctuation, *especially* in comments that will become part of the Doxygen API documentation. It is tedious for reviewers to correct your spelling -- a spell checker is your friend and ours here. We know spelling and grammatical errors will creep in, but the fewer the better. If you are not a native English speaker, please just do your best -- we'll help.

### Dates and times
The need for date and time stamps arises frequently enough, and causes enough trouble, that we want to state some general preferences here, although not specific requirements for any particular situation. Maybe this goes without saying, but just go ahead and use four digits for the year! Let’s not go through that again. Compact date stamps such as those appearing in file names and source comments should have the form yyyymmdd, e.g. 20060322 which has the distinct advantage of being sortable, with the most significant part first. Code that formats friendly dates for user consumption should avoid ambiguous formats like 7/5/2005 (July 5 in the U.S. and May 7 in Europe). Instead, use July 5, 2005 or 5 July 2005 or 2005-May-07, for example.
For binary time stamps generated programmatically, please give careful thought to time zone independence.


List of Contributors
--------------------
This is an attempt at a complete contributor list; please submit a PR or file an issue if you or someone else is missing.

Real name          | GitHub Id    | Contributions/expertise
-------------------|--------------|-------------------------
Michael Sherman    |@sherm1       |Lead developer; multibody dynamics
Peter Eastman      |@peastman     |Much early Simbody development; visualizer
Chris Dembia       |@chrisdembia  |Build, task space control, CMA optimizer, bug fixes & documentation
Thomas Uchida      |@tkuchida     |Rigid impact theory & code; documentation
Ian Stavness       |@stavness     |Computational geometry
Andreas Scholz     |@AndreasScholz|Computational geometry
José Rivero        |@j-rivero     |Build, especially for Debian
Steven Peters      |@scpeters     |Build and VectorIterator improvements
John Hsu           |@hsu          |Bug fixes; iterative solver & contact theory
Nate Koenig        |@nkoenig      |Bug fix
Ayman Habib        |@aymanhab     |Bug fixes; visualization, SWIG improvements
Ajay Seth          |@aseth1       |Mobilizer theory and code
Jack Wang          |@jmwang       |Bug fixes; visualization improvements
Tim Dorn           |@twdorn       |Bug fixes
Apoorva Rajagopal  |@apoorvar     |Xcode build fixes
Kevin Xu           |@kevunix      |Build fix
Elena Ceseracciu   |@elen4        |Improved dependency resolution
Kevin He           |@kingchurch   |Bug fixes
Paul Mitiguy       |              |Rotation class; dynamics
Matthew Millard    |@mjhmilla     |Bicubic spline
Jack Middleton     |              |Numerical methods
Christopher Bruns  |@cmbruns-hhmi |Molmodel
Randy Radmer       |              |Molmodel
Charles Schwieters |              |Provided initial multibody code 
Abhinandan Jain    |              |Underlying spatial algebra formulation
Isaac Newton       |              |F=ma, calculus, etc.


Contributor License Agreement
-----------------------------
Simbody is licensed under the very permissive [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0). Simbody users are not required to follow our noble egalitarian principles, nor to share their profits with us, nor even to acknowledge us (though they often do). When you make a contribution in any of the ways described above, you are agreeing to allow your contribution to be used under the same terms, adding no additional restrictions to the Simbody project nor requirements on Simbody users. 

Specifically, by contributing you are agreeing to the following terms:

  1. The code, text, or images you submit are your original work (you own and retain the copyright) or you otherwise have the right to submit the work.
  2. You grant the Simbody project, developers, and users a nonexclusive, irrevocable license to use your submission and any necessary intellectual property, under terms of the Apache 2.0 license.
  3. No part of your contribution is covered by a viral ("copyleft") license like GPL or LGPL.
  4. You are capable of granting these rights for the contribution.

If your contribution contains others' open source code licensed under Apache 2.0 or other non-viral license like BSD or MIT, it is probably fine. But be sure to mention that in the Pull Request you submit so we can discuss it. 
