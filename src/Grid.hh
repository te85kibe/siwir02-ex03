template< typename Type, uint Cellsize >
	class Grid
	{
	public:
		inline Grid();
		inline Grid( uint xsize, uint ysize );
	//	inline ~Grid();
	
		inline Type& operator()( uint x, uint y, uint f );
		inline Type  operator()( uint x, uint y, uint f ) const;

	private:
		uint xsize_;  // Number of nodes in x-dimension
		uint ysize_;  // Number of nodes in y-dimension
		Type* data_;  // Linearized, 1-dimensional representation
                 // of the 2D data grid
	};

// Implementation of the default constructor
template< typename Type, uint Cellsize >
Grid<Type,Cellsize>::Grid()
   : xsize_(0)
   , ysize_(0)
   , data_(0)
{}
// Implementation of the initialization constructor
template< typename Type, uint Cellsize >
Grid<Type,Cellsize>::Grid( uint xsize, uint ysize )
   : xsize_(xsize)
   , ysize_(ysize)
   , data_( new Type[Cellsize*xsize*ysize] )
{}

// Implementation of the function call operator
template< typename Type, uint Cellsize >
inline Type&
Grid<Type,Cellsize>::operator()( uint x, uint y, uint f )
{
   assert( x < xsize_ && y < ysize_ && f < Cellsize );
   return data_[y*xsize_*Cellsize+x*Cellsize+f];
}
// Implementation of the const function call operator
// ... Same as non-const version

