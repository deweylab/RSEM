// ***************************************************************************
// BamWriter.h (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_H
#define BAMWRITER_H

#include <api_global.h>
#include <BamAlignment.h>
#include <string>

namespace BamTools {

namespace Internal {
    class BamWriterPrivate;
} // namespace Internal

class API_EXPORT BamWriter {

    // constructor/destructor
    public:
        BamWriter(void);
        ~BamWriter(void);

    // public interface
    public:
        // closes the alignment archive
        void Close(void);
        // opens the alignment archive
        bool Open(const std::string& filename, 
                  const std::string& samHeader, 
                  const BamTools::RefVector& referenceSequences, 
                  bool writeUncompressed = false);
        // saves the alignment to the alignment archive
        void SaveAlignment(const BamTools::BamAlignment& al);

    // private implementation
    private:
        Internal::BamWriterPrivate* d;
};

} // namespace BamTools

#endif // BAMWRITER_H
