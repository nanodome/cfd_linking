#ifndef STREAMLINE_FILEXML_H
#define STREAMLINE_FILEXML_H

#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "XMLfile.h"
#include "streamline.h"

class streamlinefileXML: public XMLfile{

	/// Sampling time start
	double samples_start;

	/// Sampling time end
	double samples_end;

	/// Number of streamlines in the file
	int n_streams;

    /// streamlines
    std::vector<streamline>& streams;

    /// Enumerator for Streamlines XML TAGS
    enum Stream_Tags {GP, T_START, T_END, ID,
                      STREAM,
                      N_TIME_SAMPLES, TIME_SAMPLES,
                      T, P,
                      N_SPECIES, SPECIES, MOLAR_C, X, Y, Z};

    /// map for Strings and Tags
    std::map<std::string, Stream_Tags> tag_map;

    /// Tags Map Initialization
    void map_Init();

    /// Get streamline data from XML File
    /// \param XMLElement *cursor_ptr: pointer for exploring the <STREAM> sub_tree
    /// \param XMLElement *root_ptr: pointer to <STREAM> sub_tree root
    bool parse_Stream(XMLElement *cursor_ptr, XMLElement *root_ptr);

	/// Print streamline data to XML File
	/// \param XMLElement *cursor_ptr: pointer for expanding the <STREAM> sub_tree
	/// \param streamline stream: streamline to write
	void write_Stream(streamline stream, XMLElement *r_sub_tree);



public:

    /// Constructor (reads and stores the XML File)
    /// \param string path: Path of the XML file
	streamlinefileXML(std::string _path, std::vector<streamline>& _streams);

	///Return Streamlines from file
	std::vector<streamline> get_streamlines() const  { return streams; }

	/// Return start time
	double  get_start_time() const { return samples_start; }

	/// Return end sampling time
	double get_end_time() const { return samples_end; }

    /// extract Streamlines from XML file
    void read_Streamlines();

	
    /// write streamlines (defined outside the object) to file
	/// \param std::vector<streamline>& streamlines reference to streamlines container
	/// \param double _start_t: sampling start time
	/// \param double _end_t: sampling end time
	void write_Streamlines(std::vector<streamline>& _streamlines, 
						   double _start_t, double _end_t);

	/// write streamlines defined inside the object to file
	/// \param std::string _path: output file path
	void write_Streamlines(std::string _path);

    /// Print Streamlines
    void printstreamlines();

	

};

#endif // STREAMLINE_FILE

