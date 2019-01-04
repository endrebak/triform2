#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <vector>
#include <iterator>


typedef std::vector<uint32_t> intvec;
typedef std::pair <std::string, char> key;
typedef std::map<key, intvec> genome_map;


struct interval { uint32_t start; uint32_t end; };
typedef std::vector<interval> interval_vector;
typedef std::map<key, interval_vector> genome_intervals;


uint32_t compare_by_start_end_(const interval lhs, const interval rhs){
  if (lhs.start < rhs.start){
    return 1;
  } else if (rhs.start < lhs.start){
    return 0;
  } else if (lhs.end < rhs.end){
    return 1;
  } else {
    return 0;
      };

}


uint32_t start_end_equal_(const interval lhs, const interval rhs){
  if ((lhs.start == rhs.start) && (lhs.end == rhs.end)){
      return 1;
    } else {
    return 0;
  }
}




genome_map intervals_to_points(genome_intervals genome, uint32_t drop_duplicates, uint32_t read_width) {

  genome_intervals::iterator it;
  char strand;

  uint32_t i = 0;

  key chrom_strand;
  interval_vector intervals;
  genome_map genome_tags;
  uint32_t gap = 0;
  uint32_t start = 0;
  uint32_t end = 0;
  uint32_t width = 0;

  it = genome.begin();

  while(it != genome.end()){

    chrom_strand = it->first;
    strand = it->first.second;
    intervals = it->second;

    if (drop_duplicates){
      std::sort(intervals.begin(), intervals.end(), compare_by_start_end_);
      intervals.erase(unique(intervals.begin(), intervals.end(), start_end_equal_), intervals.end());
    }

    i = 0;
    for (i = 0; i < intervals.size(); i++){
      width = (intervals[i].end - (intervals[i].start + 1));
      gap = floor((read_width - width)/2);
      start = intervals[i].start - gap - 1;
      end = start + width + gap + gap;

      genome_tags[chrom_strand].push_back(start);
      genome_tags[chrom_strand].push_back(end);
    }

    genome.erase(chrom_strand);
    it++;

  }

  return genome_tags;
}



genome_map read_bed(char const* fileName, uint32_t drop_duplicates, uint32_t read_width)
{
  std::ifstream file(fileName);

  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  interval_vector intervals;
  genome_intervals genome;


  while(file >> chromosome >> start >> end >> junk >> junk >> strand)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_midpoint(genome, drop_duplicates);
}





genome_map read_bedpe(char const* fileName, uint32_t drop_duplicates, uint32_t read_width)
{
  std::ifstream file(fileName);

  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  interval_vector intervals;
  genome_intervals genome;


  while(file >> chromosome >> start >> junk >> junk >> junk >> end >> junk >> junk >> strand >> junk)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_midpoint(genome, drop_duplicates);

}

#include "gzstream.h"
genome_map read_bed_gz(char const* fileName, uint32_t drop_duplicates, uint32_t read_width)
{
  igzstream in(fileName);


  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;


  key chrom_strand;
  interval _interval;
  // intvec tags;
  interval_vector intervals;
  genome_intervals genome;

  while(in >> chromosome >> start >> end >> junk >> junk >> strand)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_midpoint(genome, drop_duplicates);

}


genome_map read_bedpe_gz(char const* fileName, uint32_t drop_duplicates, uint32_t read_width)
{

  igzstream in(fileName);

  std::string   chromosome;
  std::string   junk;
  uint32_t start;
  uint32_t end;
  char strand;

  key chrom_strand;
  interval _interval;
  interval_vector intervals;
  genome_intervals genome;

  while(in >> chromosome >> start >> junk >> junk >> junk >> end >> junk >> junk >> strand >> junk)
    {

      chrom_strand = std::make_pair(chromosome, strand);

      _interval = {start, end - 1};

      genome[chrom_strand].push_back(_interval);

    }

  return intervals_to_midpoint(genome, drop_duplicates);
}
