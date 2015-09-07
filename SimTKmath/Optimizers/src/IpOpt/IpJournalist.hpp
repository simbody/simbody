// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpJournalist.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPJOURNALIST_HPP__
#define __IPJOURNALIST_HPP__

#include "IpoptConfig.h"
#include "IpTypes.hpp"
#include "IpReferenced.hpp"
#include "IpSmartPtr.hpp"

#ifdef HAVE_CSTDARG
# include <cstdarg>
#else
# ifdef HAVE_STDARG_H
#  include <stdarg.h>
# else
#  error "don't have header file for stdarg"
# endif
#endif

#include <string>
#include <vector>

namespace Ipopt
{

  // forward declarations
  class Journal;
  class FileJournal;

  /**@name Journalist Enumerations. */
  //@{
  /** Print Level Enum. */
  enum EJournalLevel {
    J_INSUPPRESSIBLE=-1,
    J_NONE=0,
    J_ERROR,
    J_WARNING,
    J_SUMMARY,
    J_ITERSUMMARY,
    J_DETAILED,
    J_MOREDETAILED,
    J_VECTOR,
    J_MOREVECTOR,
    J_MATRIX,
    J_MOREMATRIX,
    J_ALL,
    J_LAST_LEVEL
  };

  /** Category Selection Enum. */
  enum EJournalCategory {
    J_DBG=0,
    J_STATISTICS,
    J_MAIN,
    J_INITIALIZATION,
    J_BARRIER_UPDATE,
    J_SOLVE_PD_SYSTEM,
    J_FRAC_TO_BOUND,
    J_LINEAR_ALGEBRA,
    J_LINE_SEARCH,
    J_HESSIAN_APPROXIMATION,
    J_SOLUTION,
    J_DOCUMENTATION,
    J_NLP,
    J_TIMING_STATISTICS,
    J_USER_APPLICATION /** This can be used by the user's application*/ ,
    J_LAST_CATEGORY
  };
  //@}

  /** Class responsible for all message output.
   * This class is responsible for all messaging and output.
   * The "printing" code or "author" should send ALL messages to the
   * Journalist, indicating an appropriate category and print level.
   * The journalist then decides, based on reader specified
   * acceptance criteria, which message is actually printed in which 
   * journals.
   * This allows the printing code to send everything, while the 
   * "reader" can decide what they really want to see.
   * 
   * Authors:
   * Authors use the 
   * Journals: You can add as many Journals as you like to the
   * Journalist with the AddJournal or the AddFileJournal methods. 
   * Each one represents a different printing location (or file).  
   * Then, you can call the "print" methods of the Journalist to output
   * information to each of the journals.
   * 
   * Acceptance Criteria: Each print message should be flagged 
   * appropriately with an EJournalCategory and EJournalLevel.
   * 
   * The AddFileJournal
   * method returns a pointer to the newly created Journal object
   * (if successful) so you can set Acceptance criteria for that
   * particular location.
   * 
   */
  class Journalist : public ReferencedObject
  {
  public:
    /**@name Constructor / Desructor. */
    //@{
    /** Constructor. */
    Journalist();

    /** Destructor... */
    virtual ~Journalist();
    //@}

    /**@name Author Methods.
     * These methods are used by authoring code, or code that wants
     * to report some information.
     */
    //@{
    /** Method to print a formatted string */
    void Printf(EJournalLevel level, EJournalCategory category,
                const char* format, ...) const;

    /** Method to print a long string including indentation.  The
     *  string is printed starting at the current position.  If the
     *  position (counting started at the current position) exceeds
     *  max_length, a new line is inserted, and indent_spaces many
     *  spaces are printed before the string is continued.  This is
     *  for example used during the printing of the option
     *  documentation. */
    void PrintStringOverLines(EJournalLevel level, EJournalCategory category,
                              Index indent_spaces, Index max_length,
                              const std::string& line) const;

    /** Method to print a formatted string with indentation */
    void PrintfIndented(EJournalLevel level,
                        EJournalCategory category,
                        Index indent_level,
                        const char* format, ...) const;

    /** Method to print a formatted string
     * using the va_list argument. */
    void VPrintf(EJournalLevel level,
                 EJournalCategory category,
                 const char* pformat,
                 va_list ap) const;

    /** Method to print a formatted string with indentation,
     * using the va_list argument. */
    void VPrintfIndented(EJournalLevel level,
                         EJournalCategory category,
                         Index indent_level,
                         const char* pformat,
                         va_list ap) const;

    /** Method that returns true if there is a Journal that would
     *  write output for the given JournalLevel and JournalCategory.
     *  This is useful if expensive computation would be required for
     *  a particular output.  The author code can check with this
     *  method if the computations are indeed required.
     */
    bool ProduceOutput(EJournalLevel level,
                       EJournalCategory category) const;


    /** Method that flushes the current buffer for all Journalists.
     Calling this method after one optimization run helps to avoid
     cluttering output with that produced by other parts of the
     program (e.g. written in Fortran) */
    void FlushBuffer() const;
    //@}

    /**@name Reader Methods.
     * These methods are used by the reader. The reader will setup the 
     * journalist with each output file and the acceptance
     * criteria for that file.
     *
     * Use these methods to setup the journals (files or other output).
     * These are the internal objects that keep track of the print levels 
     * for each category. Then use the internal Journal objects to
     * set specific print levels for each category (or keep defaults).
     *  
     */
    //@{
    /** Add a new journal.  The location_name is a string identifier,
     *  which can be used to obtain the pointer to the new Journal at
     *  a later point using the GetJournal method.
     *  The default_level is
     *  used to initialize the * printing level for all categories.
     */
    bool AddJournal(const SmartPtr<Journal> jrnl);

    /** Add a new FileJournal. fname is the name
     *  of the * file to which this Journal corresponds.  Use
     *  fname="stdout" * for stdout, and use fname="stderr" for
     *  stderr.  This method * returns the Journal pointer so you can
     *  set specific acceptance criteria.  It returns NULL if there
     *  was a problem creating a new Journal.    
     */
    SmartPtr<Journal> AddFileJournal(
      const std::string& location_name,    /** identifier */
      const std::string& fname,
      EJournalLevel default_level = J_WARNING
    );

    /** Get an existing journal.  You can use this method to change
     *  the acceptance criteria at runtime.
     */
    SmartPtr<Journal> GetJournal(const std::string& location_name);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    Journalist(const Journalist&);

    /** Overloaded Equals Operator */
    void operator=(const Journalist&);
    //@}

    //** Private Data Members. */
    //@{
    std::vector< SmartPtr<Journal> > journals_;
    //@}
  };

  /** Journal class (part of the Journalist implementation.). This
   *  class is the base class for all Journals. It controls the 
   *  acceptance criteria for print statements etc. Derived classes
   *  like the FileJournal - output those messages to specific locations
   */
  class Journal : public ReferencedObject
  {
  public:
    /** Constructor. */
    Journal(const std::string& name, EJournalLevel default_level);

    /** Destructor. */
    virtual ~Journal();

    /** Get the name of the Journal */
    std::string Name();

    /** Set the print level for a particular category. */
    void SetPrintLevel(
      EJournalCategory category, EJournalLevel level
    );

    /** Set the print level for all category. */
    void SetAllPrintLevels(
      EJournalLevel level
    );

    /**@name Journal Output Methods. These methods are called by the
     *  Journalist who first checks if the output print level and category
     *  are acceptable.
     *  Calling the Print methods explicitly (instead of through the 
     *  Journalist will output the message regardless of print level
     *  and category. You should use the Journalist to print & flush instead
     */
    //@{
    /** Ask if a particular print level/category is accepted by the
     * journal.
     */
    bool IsAccepted(
      EJournalCategory category, EJournalLevel level
    ) const;

    /** Print to the designated output location */
    void Print(const char* str)
    {
      PrintImpl(str);
    }

    /** Printf to the designated output location */
    void Printf(const char* pformat, va_list ap)
    {
      PrintfImpl(pformat, ap);
    }

    /** Flush output buffer.*/
    void FlushBuffer()
    {
      FlushBufferImpl();
    }
    //@}

  protected:
    /**@name Implementation version of Print methods. Derived classes
     * should overload the Impl methods.
     */
    //@{
    /** Print to the designated output location */
    virtual void PrintImpl(const char* str)=0;

    /** Printf to the designated output location */
    virtual void PrintfImpl(const char* pformat, va_list ap)=0;

    /** Flush output buffer.*/
    virtual void FlushBufferImpl()=0;
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    Journal();

    /** Copy Constructor */
    Journal(const Journal&);

    /** Overloaded Equals Operator */
    void operator=(const Journal&);
    //@}

    /** Name of the output location */
    std::string name_;

    /** vector of integers indicating the level for each category */
    Index print_levels_[J_LAST_CATEGORY];
  };


  /** FileJournal class. This is a particular Journal implementation that
   *  writes to a file for output. It can write to (stdout, stderr, or disk)
   *  by using "stdout" and "stderr" as filenames.
   */
  class FileJournal : public Journal
  {
  public:
    /** Constructor. */
    FileJournal(const std::string& name, EJournalLevel default_level);

    /** Destructor. */
    virtual ~FileJournal();

    /** Open a new file for the output location.
     *  Special Names: stdout means stdout,
     *               : stderr means stderr.
     *
     *  Return code is false only if the file with the given name
     *  could not be opened.
     */
    bool Open(const char* fname);

  protected:
    /**@name Implementation version of Print methods - Overloaded from
     * Journal base class.
     */
    //@{
    /** Print to the designated output location */
    virtual void PrintImpl(const char* str) override;

    /** Printf to the designated output location */
    virtual void PrintfImpl(const char* pformat, va_list ap) override;

    /** Flush output buffer.*/
    virtual void FlushBufferImpl() override;
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    FileJournal();

    /** Copy Constructor */
    FileJournal(const FileJournal&);

    /** Overloaded Equals Operator */
    void operator=(const FileJournal&);
    //@}

    /** FILE pointer for the output destination */
    FILE* file_;
  };
}

#endif
