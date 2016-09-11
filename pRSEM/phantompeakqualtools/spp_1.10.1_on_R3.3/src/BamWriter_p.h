// ***************************************************************************
// BamWriter_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_P_H
#define BAMWRITER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#include <BamAux.h>
#include <BGZF.h>
#include <string>
#include <vector>

namespace BamTools {
namespace Internal {

class BamWriterPrivate {

    // ctor & dtor
    public:
	BamWriterPrivate(void);
	~BamWriterPrivate(void);

    // "public" interface to BamWriter
    public:
	void Close(void);
	bool Open(const std::string& filename,
		  const std::string& samHeader,
		  const BamTools::RefVector& referenceSequences,
		  bool isWriteUncompressed);
	void SaveAlignment(const BamAlignment& al);

    // internal methods
    public:
	const unsigned int CalculateMinimumBin(const int begin, int end) const;
	void CreatePackedCigar(const std::vector<BamTools::CigarOp>& cigarOperations, std::string& packedCigar);
	void EncodeQuerySequence(const std::string& query, std::string& encodedQuery);

    // data members
    public:
	BgzfData mBGZF;
	bool IsBigEndian;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMWRITER_P_H
