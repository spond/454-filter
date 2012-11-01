
#define DEFAULT_SZ 1024L

// vec_t ------------------------------------------------------------------------------------------------------------ //

template<class T>
inline
void vec_t<T>::__resize( long rlen )
{
    // don't forget the null terminating byte
    rlen += 1;

    if ( rlen > capacity ) {
        const long incr = ( capacity / 8 < DEFAULT_SZ ) ? DEFAULT_SZ : capacity / 8,
                   rcap = incr * ( rlen / incr + 1 );
        void * ptr = realloc( data, sizeof( T ) * rcap );
        __CHECK_PTR( ptr )
        data = reinterpret_cast<T *>( ptr );
        capacity = rcap;
    }
}

template<class T>
void vec_t<T>::init()
{
    void * ptr = NULL;

    if ( !data )
        ptr = malloc( sizeof( T ) * DEFAULT_SZ );
    else if ( capacity > DEFAULT_SZ )
        ptr = realloc( data, sizeof( T ) * DEFAULT_SZ );
    else
        goto end;

    __CHECK_PTR( ptr )
    data = reinterpret_cast<T *>( ptr );
    capacity = DEFAULT_SZ;
end:
    len = 0;
    data[len] = T( NULL );
}

template<class T>
vec_t<T>::vec_t() : data( NULL )
{
    init();
}

template<class T>
vec_t<T>::~vec_t()
{
    free( data );
}

template<class T>
void vec_t<T>::append( const T item )
{
    __resize( len + 1 );
    data[len] = item;
    len += 1;

    if ( len >= capacity )
        fprintf( stderr, "len >= capacity" );

    data[len] = T( NULL );
}

template<class T>
void vec_t<T>::clear()
{
    // init resizes capacity to DEFAULT_SZ
    init();
}

template<class T>
int vec_t<T>::compare( vec_t<T> & vec ) const
{
    long minlen, i;

    if ( len <= vec.length() )
        minlen = len;
    else
        minlen = vec.length();

    for ( i = 0; i < minlen; ++i ) {
        if ( data[i] < vec[i] )
            return -1;
        else if ( data[i] > vec[i] )
            return 1;
    }

    if ( len == vec.length() )
        return 0;

    return ( len < vec.length ) ? -1 : 1;
}

template<class T>
void vec_t<T>::extend( const T * vec, long nitem )
{
    if ( nitem <= 0 )
        return;
    else {
        __resize( len + nitem );
        memcpy( data + len, vec, sizeof( T ) * nitem );
        len += nitem;
        data[len] = T( NULL );
    }
}

template<class T>
void vec_t<T>::extend( vec_t<T> & vec, long from, long to )
{
    long nitem = to - from;

    if ( vec.length() <= 0 || nitem <= 0 )
        return;
    else
        extend( vec.data + from, nitem );
}

template<class T>
void vec_t<T>::extend( vec_t<T> & vec )
{
    extend( vec.data, vec.length() );
}

template<class T>
long vec_t<T>::length() const
{
    return len;
}

template<class T>
void vec_t<T>::sort()
{
    qsort( data, len, sizeof( T ), __elem_cmp<T> );
}

template<class T>
T vec_t<T>::operator[]( const long i ) const
{
    return data[i];
}
