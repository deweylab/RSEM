// ***************************************************************************
// BamReader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 November 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#ifndef BAMREADER_P_H
#define BAMREADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include <BamAlignment.h>
#include <BamIndex.h>
#include <BGZF.h>
#include <string>

namespace BamTools {

class BamReader;

namespace Internal {

class BamReaderPrivate {

    // enums
    public: enum RegionState { BEFORE_REGION = 0
			     , WITHIN_REGION
			     , AFTER_REGION
			     };

    // ctor & dtor
    public:
	BamReaderPrivate(BamReader* parent);
	~BamReaderPrivate(void);

    // 'public' interface to BamReader
    public:

	// file operations
	void Close(void);
	bool Open(const std::string& filename,
		  const std::string& indexFilename,
		  const bool lookForIndex,
		  const bool preferStandardIndex);
	bool Rewind(void);
	bool SetRegion(const BamRegion& region);

	// access alignment data
	bool GetNextAlignment(BamAlignment& bAlignment);
	bool GetNextAlignmentCore(BamAlignment& bAlignment);

	// access auxiliary data
	const std::string GetHeaderText(void) const;
	int GetReferenceID(const std::string& refName) const;

	// index operations
	bool CreateIndex(bool useStandardIndex);
	void SetIndexCacheMode(const BamIndex::BamIndexCacheMode mode);

    // 'internal' methods
    public:

	// ---------------------------------------
	// reading alignments and auxiliary data

	// adjusts requested region if necessary (depending on where data actually begins)
	void AdjustRegion(BamRegion& region);
	// fills out character data for BamAlignment data
	bool BuildCharData(BamAlignment& bAlignment);
	// checks to see if alignment overlaps current region
	RegionState IsOverlap(BamAlignment& bAlignment);
	// retrieves header text from BAM file
	void LoadHeaderData(void);
	// retrieves BAM alignment under file pointer
	bool LoadNextAlignment(BamAlignment& bAlignment);
	// builds reference data structure from BAM file
	void LoadReferenceData(void);
	// mark references with 'HasAlignments' status
	void MarkReferences(void);

	// ---------------------------------
	// index file handling

	// clear out inernal index data structure
	void ClearIndex(void);
	// loads index from BAM index file
	bool LoadIndex(const bool lookForIndex, const bool preferStandardIndex);

    // data members
    public:

	// general file data
	BgzfData  mBGZF;
	std::string HeaderText;
	BamIndex* Index;
	RefVector References;
	bool      HasIndex;
	int64_t   AlignmentsBeginOffset;
	std::string    Filename;
	std::string    IndexFilename;

//	Internal::BamHeader* m_header;

	// index caching mode
	BamIndex::BamIndexCacheMode IndexCacheMode;

	// system data
	bool IsBigEndian;

	// user-specified region values
	BamRegion Region;
	bool HasAlignmentsInRegion;

	// parent BamReader
	BamReader* Parent;

	// BAM character constants
	const char* DNA_LOOKUP;
	const char* CIGAR_LOOKUP;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMREADER_P_H
