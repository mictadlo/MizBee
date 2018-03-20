import processing.opengl.*;
import processing.pdf.*;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.util.Collections;

Genome _genome;
GenomeView _genome_viewer;
ChromosomeView _chromosome_viewer;
BlockView _block_viewer;
PFont sans_15, sans_12, sans_small, sans_large;

float _radius;

float BACK_COLOR = 224;
float FILL_COLOR = 165;
color _background = color( BACK_COLOR );
//color _background = color( 255 );
color _fill_color = color( FILL_COLOR );
int CHROM_COLOR = 165;

int _w, _h;
int _view_spacing;
int[] _genome_view_origin, _genome_view_x_partition;
int[] _chromosome_view_origin, _chromosome_view_x_partition;
int[] _block_view_origin, _block_view_x_partition;

int _genome_viewer_w, _chromosome_viewer_w, _chromosome_viewer_h, _block_viewer_w, _block_viewer_h;

boolean _over_start_shutter = false;
boolean _over_end_shutter = false;
boolean _last_shutter_start = false;
boolean _over_chrom_start_shutter = false;
boolean _over_chrom_end_shutter = false;
boolean _last_chrom_shutter_start = false;

boolean _scrolling = false;

boolean _print_pdf = false;
int _pdf_counter = 0;

void setup() 
{ 
  // set the window size
  int offset = 50;

  _w = displayWidth;
  _h = displayHeight;

  //_w = 1350;
  //_h = 900;

  _w = min( 1200, _w);
  _h = min( 800, _h );
//  _w = 2000;
//  _h = 1100;

  //_w = 800;
  //_h = 600;

  //_w = 1900;
  //_h = 1200;
  

  _w -= offset;
  _h -= offset;  

  size( _w, _h, OPENGL );

  _genome = new Genome( "./config.txt" );


  _genome_view_origin = new int[2];
  _genome_view_x_partition = new int[2];
  _chromosome_view_origin = new int[2];
  _chromosome_view_x_partition = new int[2];
  _block_view_origin = new int[2];
  _block_view_x_partition = new int[2];

  setDimensions( _w, _h );

  _genome_viewer = new GenomeView( _genome, _genome_viewer_w );
  _chromosome_viewer = new ChromosomeView( _genome, _chromosome_viewer_w, _chromosome_viewer_h );
  _block_viewer = new BlockView( _genome, _block_viewer_w, _block_viewer_h );


  sans_15 = createFont("Verdana", 14, true); 
  sans_12 = createFont("Verdana", 12, true); 
  sans_small = createFont( "Verdana", 12, true);
  sans_large = createFont( "Verdana", 18, true);

  //  frame.setResizable(true);
  //  frame.addComponentListener(new ComponentAdapter() {
  //   public void componentResized(ComponentEvent e) {
  //     if(e.getSource()==frame) 
  //     { 
  //       _w = frame.getWidth();
  //       _h = frame.getHeight();
  //       frame.setSize(_w,_h); 
  //       _h -= 22; // offset on mac
  //       setDimensions( _w, _h );
  //       setViewerDimensions();
  //     }
  //   }
  //  } );
} 

void setDimensions( int w, int h ) {
  _w = w;
  _h = h;

  _view_spacing = int((float)_w*0.05);

  _block_viewer_w = int((float)_w*0.1);
  _block_viewer_h = int((float)_h*0.75);

  _chromosome_viewer_w = int((float)_w*0.125);
  _chromosome_viewer_h = int((float)_h*0.85);

  int remaining_space = _w - _block_viewer_w - _chromosome_viewer_w - int(2.5*_view_spacing) - _view_spacing/3;

  //_genome_viewer_w = int((float)_w/1.75);
  _genome_viewer_w = min( remaining_space, _h );
  //_genome_viewer_w =1060;

  _block_viewer_w = _w - _genome_viewer_w - _chromosome_viewer_w - int(2.5*_view_spacing) - _view_spacing/3;

  _genome_view_origin[0] = _view_spacing/6;
  //_genome_view_origin[0] = (_w-_genome_viewer_w)/2;
  _genome_view_origin[1] = int(0.5*(_h - _genome_viewer_w));
  _genome_view_x_partition[0] = _genome_view_origin[0];
  _genome_view_x_partition[1] = _genome_view_x_partition[0] + _genome_viewer_w;


  _chromosome_view_origin[0] = _genome_view_origin[0] + _genome_viewer_w + _view_spacing;
  _chromosome_view_origin[1] = int(0.5*(_h - _chromosome_viewer_h));
  _chromosome_view_x_partition[0] = _chromosome_view_origin[0];
  _chromosome_view_x_partition[1] = _chromosome_view_x_partition[0] + _chromosome_viewer_w;


  _block_view_origin[0] = _chromosome_view_origin[0] + _chromosome_viewer_w + _view_spacing;
  _block_view_origin[1] = int(0.45*(_h - _block_viewer_h));
  _block_view_x_partition[0] = _block_view_origin[0];
  _block_view_x_partition[1] = _block_view_x_partition[0] + _block_viewer_w;
}

void setViewerDimensions()
{
  _genome_viewer.setDimensions( _genome_viewer_w );
  _chromosome_viewer.setDimensions( _chromosome_viewer_w, _chromosome_viewer_h );
  _block_viewer.setDimensions( _block_viewer_w, _block_viewer_h );
}

void draw() { 

  if ( _print_pdf ) {
    String filename = "output_" + nf(_pdf_counter++, 2) + ".pdf";
    beginRecord(PDF, filename);
  }

  boolean mizbee = true;

  if ( mizbee ) {
    background( _background ); 

    smooth();

    pushMatrix();
    translate( _genome_view_origin[0], _genome_view_origin[1] );
    _genome_viewer.render();
    popMatrix();

    pushMatrix();
    translate( _chromosome_view_origin[0], _chromosome_view_origin[1] );
    _chromosome_viewer.render();
    popMatrix();

    pushMatrix();
    translate( _block_view_origin[0], _block_view_origin[1] );
    _block_viewer.render();
    popMatrix();
  }

  else {


    //noSmooth();
    background( color(255) );

    pushMatrix();
    translate( 150, 50 );

    long pw = 600;
    long ph = 600;

    long max_vx = 300000000;
    long max_vy = 300000000;
    
    _selected_chromosome = 0;
    
    long steps_x = 7;
    long steps_y = 7;
    long delta_x = max_vx / (steps_x-1);
    long delta_y = max_vy / (steps_y-1);
    
    long pix_delta_x = pw / (steps_x - 1);
    long pix_delta_y = ph / (steps_y - 1);

    color line_col = color(50);
    stroke( line_col );
    noFill();
    strokeWeight( 2 );
    rect( 0, 0, pw, ph );
    
    textFont( sans_small );
    stroke( line_col );
    fill( line_col );
    
    long location = 0;
    long len = 10;
    textAlign( CENTER, TOP );
    for ( int i = 0; i < steps_x; i++ ) {
      line( location, ph, location, ph+len );
      text( nf((i*delta_x), 1,0), location, ph+len*1.5 );
      
      location += pix_delta_x;
    }
    
    location = 0;
    textAlign( RIGHT, CENTER );
    for ( int i = 0; i < steps_x; i++ ) {
      line( 0, ph-location, -len, ph-location );
      text( nf((i*delta_x), 1,0), -len*1.5, ph-location );
      
      location += pix_delta_y;
    }
    
    textFont( sans_large );
    textAlign( CENTER, CENTER );
    text( "Human chrY", pw/2, ph + len*5 );
    
    pushMatrix();
    rotate( -PI/2.0 );
    text( "Lizard chr1", -pw/2, -len*11 );
    popMatrix();
    
    float pos_x, pos_y;
    for ( int i = 0; i < _genome.getChromosome(_selected_chromosome).numSubunits(); i++ ) {
      pos_x = _genome.getChromosome(_selected_chromosome).getBlock(i)._start;
      pos_y = _genome.getChromosome(_selected_chromosome).getBlock(i).pair()._start;
      
      pos_x = ((float)pos_x/(float)delta_x)*(float)pix_delta_x;
      pos_y = ((float)pos_y/(float)delta_y)*(float)pix_delta_y;
      
      strokeWeight( 4 );
      point( pos_x, pos_y );
    }

    popMatrix();
  }

  if ( _print_pdf ) {
    _print_pdf = false;
    endRecord();
  }
} 

void mousePressed() {
  if ( (mouseX >= _genome_view_x_partition[0]) && (mouseX <= _genome_view_x_partition[1]) ) {
    float x = (float)mouseX - (float)_genome_view_origin[0] - _genome_viewer.centerX();
    float y = (float)mouseY - (float)_genome_view_origin[1] - _genome_viewer.centerY();

    _over_start_shutter = _genome_viewer.overStartShutter( x, y );
    _over_end_shutter = _genome_viewer.overEndShutter( x, y );

    if ( _over_start_shutter && _over_end_shutter ) {
      _over_start_shutter = _last_shutter_start;
      _over_end_shutter = !_last_shutter_start;
    }
    else if ( _over_start_shutter || _over_end_shutter ) {
      _last_shutter_start = _over_start_shutter;
    }
    else {
      if ( _genome_viewer.overChromosomeNumber( x, y ) ) {
        _chromosome_viewer.setScale();
        float start = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowStartDegree());
        float end = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowEndDegree());
        _chromosome_viewer.setWindowStartAndEnd( start, end );
      }
      else {
        if ( _genome_viewer.blockSelected( x, y ) ) {
          _block_viewer.setUpZoom(_genome.getChromosome(_selected_chromosome).getBlock(_selected_block).minGeneBlockRatio());
        }
        _block_viewer.setGeneSelected( false );
      }
    }
  } 
  else if ( (mouseX >= _chromosome_view_x_partition[0]) && (mouseX <= _chromosome_view_x_partition[1]) ) {
    float x = (float)mouseX - (float)_chromosome_view_origin[0];
    float y = (float)mouseY - (float)_chromosome_view_origin[1];

    if ( _chromosome_viewer.overWindowShutterArea((int)x, (int)y) ) {
      _over_chrom_start_shutter = _chromosome_viewer.overStartShutter( x, y );
      _over_chrom_end_shutter = _chromosome_viewer.overEndShutter( x, y );

      if ( _over_chrom_start_shutter && _over_chrom_end_shutter ) {
        _over_chrom_start_shutter = _last_chrom_shutter_start;
        _over_chrom_end_shutter = !_last_chrom_shutter_start;
      }
      else if ( _over_chrom_start_shutter || _over_chrom_end_shutter ) {
        _last_chrom_shutter_start = _over_chrom_start_shutter;
      }
    }
    else if ( _chromosome_viewer.overTextBox( (int)x, (int)y ) )
    {
    }
    else {
      _over_chrom_start_shutter = _over_chrom_end_shutter = false;

      if ( _chromosome_viewer.blockSelected( (mouseX-_chromosome_view_origin[0]), (mouseY-_chromosome_view_origin[1]) ) ) {
        _block_viewer.setUpZoom(_genome.getChromosome(_selected_chromosome).getBlock(_selected_block).minGeneBlockRatio());
      }
      _block_viewer.setGeneSelected( false );
    }
  }
  else if ( (mouseX >= _block_view_x_partition[0]) /*&& (mouseX <= _block_view_x_partition[1])*/ ) {
    int x = (mouseX-_block_view_origin[0]);
    int y = (mouseY-_block_view_origin[1]);
    if ( !_block_viewer.invertSelected( x, y ) ) {
      if ( !_block_viewer.overZoomControl(x, y) ) {
        _scrolling = _block_viewer.overScrollControl( x, y );
        if ( !_scrolling ) {
          if (mouseEvent.getClickCount()==2) {
            _block_viewer.scroll( y );
          }
          else _block_viewer.geneSelected( x, y );
        }
      }
    }
  }
  else {
    _block_viewer.setGeneSelected( false );
  }
}

void mouseReleased() {
  _scrolling = false;
}

void mouseDragged() {

  if ( (mouseX >= _genome_view_x_partition[0]) && (mouseX <= _genome_view_x_partition[1]) ) {
    float x = (float)mouseX - (float)_genome_view_origin[0] - _genome_viewer.centerX();
    float y = (float)mouseY - (float)_genome_view_origin[1] - _genome_viewer.centerY();
    if ( _over_start_shutter || _over_end_shutter ) {
      _genome_viewer.updateShutter( x, y, _over_start_shutter, _over_end_shutter, _genome.getChromosome(_selected_chromosome).startDegree(), 
      _genome.getChromosome(_selected_chromosome).stopDegree() );

      float start = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowStartDegree());
      float end = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowEndDegree());
      _chromosome_viewer.setWindowStartAndEnd( start, end );
    }
    else {
      _genome_viewer.overAlphaSlider( (int)x, (int)y );
    }
  }
  else if ( (mouseX >= _chromosome_view_x_partition[0]) && (mouseX <= _chromosome_view_x_partition[1]) ) {
    if ( _over_chrom_start_shutter || _over_chrom_end_shutter ) {
      float x = (float)mouseX - (float)_chromosome_view_origin[0];
      float y = (float)mouseY - (float)_chromosome_view_origin[1];

      _chromosome_viewer.updateShutter( x, y, _over_chrom_start_shutter, _over_chrom_end_shutter );

      float start = _genome.getChromosome(_selected_chromosome).convertRelativePositionToDegree(_chromosome_viewer.windowStart());
      float end = _genome.getChromosome(_selected_chromosome).convertRelativePositionToDegree(_chromosome_viewer.windowEnd());
      _genome_viewer.setWindowStartAndEnd( start, end );
    }
  }
  else if ( (mouseX >= _block_view_x_partition[0]) /*&& (mouseX <= _block_view_x_partition[1])*/ ) {
    int x = (mouseX-_block_view_origin[0]);
    int y = (mouseY-_block_view_origin[1]);
    if ( _scrolling ) {
      _block_viewer.scrollTheBar( y );
    }
    else {
      _block_viewer.overZoomControl(x, y);
    }
  }
}

void mouseMoved() {

  if ( (mouseX >= _block_view_x_partition[0]) && (mouseX <= _block_view_x_partition[1]) ) {
    int x = mouseX - _block_view_origin[0];
    int y = mouseY - _block_view_origin[1];

    _block_viewer.overGene( x, y );
  }
  else {
    _gene_rolled_over = false;
  }

  if ( (mouseX >= _chromosome_view_x_partition[0]) && (mouseX <= _chromosome_view_x_partition[1]) ) {
    int x = mouseX - _chromosome_view_origin[0];
    int y = mouseY - _chromosome_view_origin[1];

    _chromosome_viewer.overBlock( x, y );
  }
}

void keyPressed() {
  println( "HERE" );
  if ( key == CODED ) {
    if ( keyCode == UP ) {
      _chromosome_viewer.upKeyPressed();
      _block_viewer.setUpZoom(_genome.getChromosome(_selected_chromosome).getBlock(_selected_block).minGeneBlockRatio());
    } 
    else if (keyCode == DOWN) {
      _chromosome_viewer.downKeyPressed();
      _block_viewer.setUpZoom(_genome.getChromosome(_selected_chromosome).getBlock(_selected_block).minGeneBlockRatio());
    } 
    else if (keyCode == LEFT) {
      _genome_viewer.leftKeyPressed();
      _chromosome_viewer.setScale();

      float start = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowStartDegree());
      float end = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowEndDegree());
      _chromosome_viewer.setWindowStartAndEnd( start, end );
    } 
    else if (keyCode == RIGHT) {
      _genome_viewer.rightKeyPressed();
      _chromosome_viewer.setScale();

      float start = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowStartDegree());
      float end = _genome.getChromosome(_selected_chromosome).convertDegreeToRelativePosition(_genome_viewer.windowEndDegree());
      _chromosome_viewer.setWindowStartAndEnd( start, end );
    }
  } 
  else if ( (key >= 48) && (key <= 57 ) ) {
    _chromosome_viewer.coordInput( key );
  }
  else if (key == BACKSPACE ) {
    _chromosome_viewer.removeOneText();
  }
  else if ( (key == RETURN) || (key == ENTER) ) {
    _chromosome_viewer.returnPressed();
    _block_viewer.setUpZoom(_genome.getChromosome(_selected_chromosome).getBlock(_selected_block).minGeneBlockRatio());
  }
  else if ( key == 'P' ) {
    _print_pdf = true;
  }
}





