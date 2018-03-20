int _selected_genome = 0;
int MAX_NUM_TRACKS = 4;

boolean _have_tracks = false;

class Genome {
  Chromosome[] _chromosome;
  int[] _num_chromosomes;
  long[] _length;
  String[] _genome_tag;
  int[] _genome_chrom_start_index;
  boolean _two_genomes;
  float _min_gene_block_ratio;
  int _num_tracks;
  String[] _track_tags;

  Genome( String config_file ) {
    String[] rows = loadStrings(config_file);
      
    buildChromosomesGenesAndBlocks( rows[0] );
    
    if ( rows.length > 2 ) {
      _have_tracks = true;
      
      // get the number of tracks
      String[] cols = splitTokens(rows[1]);
      _num_tracks = min(parseInt( cols[1] ), MAX_NUM_TRACKS); // cap the number of tracks at 4 for now
      _track_tags = new String[_num_tracks];
      
      for ( int i = 0; i < _num_chromosomes[0]; i++ ) _chromosome[i].setNumTracks( _num_tracks );
      
      for ( int i = 0; i < _num_tracks; i++ ) {
        cols = splitTokens( rows[i+2] );
        
        readInTrack( cols[0], i );
        _track_tags[i] = cols[1];
      }
    }
  }
  
  int numTracks() {
    return _num_tracks;
  }
  
  String trackTag( int t ) {
    return _track_tags[t];
  }

  long len( int which_genome ) {
    return _length[which_genome];
  }

  int numChromosomes( int which_genome ) {
    return _num_chromosomes[which_genome];
  }
  
  int genomeChromStartIndex( int which_genome) {
    return _genome_chrom_start_index[which_genome];
  }

  Chromosome getChromosome( int index ) {
    return _chromosome[index];
  }
  
  float minGeneBlockRatio() {
    return _min_gene_block_ratio;
  }

  int findChromosomeIndex( String tag ) {
    if ( _two_genomes ) {
      for ( int i = 0; i < _chromosome.length; i++ ) {
        if ( tag.equals(_chromosome[i].tag()) ) {
          return i;
        }  
      }
    }
    else {
      for ( int i = _genome_chrom_start_index[1]; i < _chromosome.length; i++ ) {
        if ( tag.equals(_chromosome[i].tag()) ) {
          return i;
        }  
      }
    }
    println( "COULDN'T FIND CHROMOSOME TAG: " + tag );
    return -1;
  }
  
  int findChromosomeIndex( String tag, boolean search_second_half ) {
    if ( _two_genomes && !search_second_half ) {
      for ( int i = 0; i < _chromosome.length; i++ ) {
        if ( tag.equals(_chromosome[i].tag()) ) {
          return i;
        }  
      }
    }
    else {
      for ( int i = _genome_chrom_start_index[1]; i < _chromosome.length; i++ ) {
        if ( tag.equals(_chromosome[i].tag()) ) {
          return i;
        }  
      }
    }
    println( "COULDN'T FIND CHROMOSOME TAG: " + tag );
    return -1;
  }

  boolean parseOrientation( String s ) {
    if ( (s.substring(0,1)).equals("+") ) {
      s = s.substring(1, s.length());
    }
    if( parseInt(s) == 1 ) {
      return true;
    }
    else {
      return false;
    } 
  }
  
  void readInTrack( String file, int index ) {
    String[] rows = loadStrings( file );
    
    String[] cols;
    Annotation a;
    int c;
    for ( int i = 0; i < rows.length; i++ ) {
      cols = splitTokens( rows[i], ",\t" );
      
      if ( cols[0].equals( "Un" ) ) continue;
      
      // FOR SINGLE GE
      c = findChromosomeIndex( cols[0] );
      if ( !_two_genomes ) c -= _num_chromosomes[0];
      
      a = new Annotation( parseInt(cols[1]), parseInt(cols[2]), _chromosome[c].len(), c );
      a.setTag( cols[3] );
      
      _chromosome[c].addAnnotation( index, a );
    }
  }

  void buildChromosomesGenesAndBlocks( String file ){
    String[] rows = loadStrings(file);

    // get the tag of the first genome
    String[] two_cols = splitTokens(rows[0], ",\t");
    if ( !two_cols[0].equals("genome") ) {
      println( "FILE FORMAT ERROR: expecting 'genome' " + rows[0] + " " + two_cols[0] );
      exit();
    }

    _genome_tag = new String[2];
    _num_chromosomes = new int[2];
    _genome_chrom_start_index = new int[2];
    _length = new long[2];

    _genome_tag[0] = two_cols[1];

    // get the number of chromosomes
    two_cols = splitTokens(rows[1]);
    if ( !two_cols[0].equals( "chromosomes" ) ) {
      println( "FILE FORMAT ERROR: expecting 'chromosomes' " );
      exit();
    }

    _num_chromosomes[0] = parseInt( two_cols[1] );
    _genome_chrom_start_index[0] = 0;

    // get the tag of the second genome (if it exists!)
    two_cols = splitTokens(rows[3+_num_chromosomes[0]], ",\t");
    if ( two_cols[0].equals("genome") ) {
      _genome_tag[1] = two_cols[1];

      // get the number of chromosomes
      two_cols = splitTokens(rows[4+_num_chromosomes[0]]);
      if ( !two_cols[0].equals( "chromosomes" ) ) {
        println( "FILE FORMAT ERROR: expecting 'chromosomes' " );
        exit();
      }
      _num_chromosomes[1] = parseInt( two_cols[1] );

      _two_genomes = true;
    }
    else {
      _genome_tag[1] = _genome_tag[0];
      _num_chromosomes[1] = _num_chromosomes[0];

      _two_genomes = false;
    }

    _genome_chrom_start_index[1] = _num_chromosomes[0];
    _chromosome = new Chromosome[_num_chromosomes[0]+_num_chromosomes[1]];

    println( _genome_tag[0] + ": " + _num_chromosomes[0] + " chromosomes");
    if ( _two_genomes ) 
      println( _genome_tag[1] + ": " + _num_chromosomes[1] + " chromosomes");

    // get the chromosome lengths    
    long[] l = new long[_num_chromosomes[0]];
    String[] t = new String[_num_chromosomes[0]];
    int counter = 1;
    Integer i_val;
    for ( int i = 0; i < _num_chromosomes[0]; i++ ) {
      two_cols = splitTokens(rows[++counter]);
      t[i] = two_cols[0];
      i_val = parseInt(two_cols[1]);
      l[i] = i_val.longValue();
    }

    // get the total length of the genome
    _length[0] = 0;
    for ( int i = 0; i < _num_chromosomes[0]; i++ ) {
      _length[0] += l[i];
    }

    // and create the chromosomes for the first genome
    int start = 1;
    for ( int i = 0; i < _num_chromosomes[0]; i++ ) {
      _chromosome[i] = new Chromosome( start, start+l[i], _length[0], i );
      _chromosome[i].setTag( t[i] );
     
      start += l[i];
    }

    if ( _two_genomes ) {      
      l = new long[_num_chromosomes[1]];
      t = new String[_num_chromosomes[1]];
      counter = 4+_num_chromosomes[0];
      for ( int i = 0; i < _num_chromosomes[1]; i++ ) {
        two_cols = splitTokens(rows[++counter]);
        t[i] = two_cols[0];
        l[i] = parseInt(two_cols[1]);      
      }

      // get the total length of the genome
      _length[1] = 0;
      for ( int i = 0; i < _num_chromosomes[1]; i++ ) {
        _length[1] += l[i];
      }

      // and create the chromosomes for the first genome
      start = 1;
      for ( int i = 0; i < _num_chromosomes[1]; i++ ) {
        _chromosome[i+_genome_chrom_start_index[1]] = new Chromosome( start, start+l[i], _length[1], i+_genome_chrom_start_index[1] );
        _chromosome[i+_genome_chrom_start_index[1]].setTag( t[i] );

        start += l[i];
      }
    }
    else {
      _length[1] = _length[0];
      for ( int i = 0; i < _num_chromosomes[1]; i++ ) {
        _chromosome[i+_genome_chrom_start_index[1]] = new Chromosome( _chromosome[i].start(), _chromosome[i].stop(), _length[1], i+_genome_chrom_start_index[1] );
        _chromosome[i+_genome_chrom_start_index[1]].setTag( _chromosome[i].tag() );
      }
    }

    // determine the number of genes 
    int num_genes = 0;
    int num_blocks = 0;
    int c1, c2, ng;
    Block b1, b2;
    boolean o;
    String orie3nt;
    String[] cols;

    int i = 3 + _num_chromosomes[0];
    if ( _two_genomes ) {
      i += 3 + _num_chromosomes[1];
    }
    
    ArrayList blocks = new ArrayList();

    _min_gene_block_ratio = Float.MAX_VALUE;
    
    String s1, s2;
    while ( i < rows.length ) {
      cols = splitTokens(rows[i], ",\t");
      
      // get the two chromosomes
      s1 = cols[0];
      s2 = cols[4];
      
      if ( !_two_genomes ) {
        c1 = findChromosomeIndex( s1 ) - _num_chromosomes[0];
      }
      else c1 = findChromosomeIndex( s1 );
      c2 = findChromosomeIndex( s2, true );
      ng = parseInt( cols[9] );

      // first, build the block and its pair
      b1 = new Block( parseInt(cols[1]), parseInt(cols[2]), _chromosome[c1].len(), c1 );
      b1.setOrientation( parseOrientation(cols[3]) );
      b2 = new Block( parseInt(cols[5]), parseInt(cols[6]), _chromosome[c2].len(), c2 );    
      b2.setOrientation( parseOrientation(cols[7]) );

      b1.setPairOrientation( b1.orientation() == b2.orientation() );
      b1.setPairValue( parseFloat(cols[8]) );
      
      // now, add the genes to both
      Gene g1, g2;    
      for ( int j = 1; j <= ng; j++ ) {
        cols = splitTokens(rows[i+j], ",\t");

        g1 = new Gene( parseInt(cols[0]), parseInt(cols[1]), _chromosome[c1].len(), c1 );
        g1.setOrientation( parseOrientation(cols[2]) );   
        if ( cols.length > 7 ) 
          g1.setTag( cols[7] );     

        g2 = new Gene( parseInt(cols[3]), parseInt(cols[4]), _chromosome[c2].len(), c2 );
        g2.setOrientation( parseOrientation(cols[5]) ); 
       if ( cols.length > 7 )   
          g2.setTag( cols[8] );     

        g1.setPair( g2 );
        g1.setPairValue( parseFloat(cols[6]) );
        g1.setPairOrientation( (g1.orientation() == g2.orientation()) );

        b1.addGene( g1 );
        b2.addGene( g2 );
        
        // find the min gene/block length ratio
        _min_gene_block_ratio = min( min(_min_gene_block_ratio, (float)g1.len()/(float)b1.len()), (float)g2.len()/(float)b2.len() );
      }
      num_genes += ng; 

      b1.setPair( b2 );
      
      // find the min gene block length ratio
      b1._min_gene_block_ratio = Float.MAX_VALUE;
      for ( int g = 0; g < b1.numSubunits(); g++ ) {
        b1._min_gene_block_ratio = min( min( b1._min_gene_block_ratio, (float)b1.getGene(g).len()/(float)b1.len() ), (float)b1.getGene(g).pair().len()/(float)b2.len() );
      }

      // add the block to the appropriate chromosome
      //_chromosome[c1].addBlock( b1 );
      blocks.add( b1 );
     
      i += ng+2; 
      ++num_blocks;      
    } 
       
    Collections.sort( blocks );    
    Block bbb;
    for ( int j = 0; j < num_blocks; j++ ) {
      bbb = (Block)blocks.get(j);
      _chromosome[bbb.refNumber()].addBlock(bbb);
    }
    
    println( "number of blocks = " + num_blocks );
    println( "number of genes = " + num_genes );
  }
  


}










