from future.utils import iteritems
import sqlite3, util, os, cPickle as pkl, dateutil, pandas as pd

def create_table(conn, create_table_sql):
    """ create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.executescript(create_table_sql)
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')
        
def table_exists(conn, table_name):
    """
    Check if table exists.
    Parameters
    ----------
    conn : db connection
    table_name : str
        Name of table to drop

    Returns
    -------
    bool
        True if table exists.
    """
    try:
        c = conn.cursor()
        c.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{}'".format(table_name))
        #if the count is 1, then table exists
        if c.fetchone()[0]==1:
            return True
        else:
            return False
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')
        raise e

def drop_table(conn, table_name):
    """
    Drop a certain table name.

    Parameters
    ----------
    conn : db connection
    table_name : str
        Name of table to drop.
    """
    try:
        c = conn.cursor()
        c.execute("SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{}'".format(table_name))
        #if the count is 1, then table exists and should be dropped
        if c.fetchone()[0]==1:
            c.executescript("drop table {};".format(table_name))
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')

def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    if not os.path.exists(db_file):
        util.clrd_print("DB '%s' does not exist. Creating new DB."%db_file, "warn")

    conn = None
    try:
        conn = sqlite3.connect(db_file)
        # SELECT returns a list of dict with column values
        conn.row_factory = sqlite3.Row
        # enforce foreign keys
        c = conn.cursor()
        c.execute("PRAGMA foreign_keys=ON;")
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')
 
    return conn

def add_record(conn, table, mode, record):
    """
    Add a record to the DB.
    Parameters
    ----------
    conn : db connection

    table : str
        Name of DB table.
    
    record : dict
        Dict with keys matching table column names and value types matching column types.

    mode : str
        Insertion mode. Choose from 'insert', 'insert-or-ignore'.
        
    Returns
    -------
    int
        If new record is inserted, returns the record (row) DB id, otherwise, return existing record ID.
    """
    try:
        c = conn.cursor()
        if mode == 'insert':
            sql_mode = 'INSERT'
        elif mode == 'insert-or-ignore':
            sql_mode = 'INSERT OR IGNORE'
        else:
            raise ValueError("Invalid record addition mode '{}'.\n".format(mode))

        # do not insert NULL values
        filtered_record = {key:val for key, val in iteritems(record) if val is not None}
        sql_query = '%s INTO %s (%s) VALUES (%s)' % (sql_mode, table, ', '.join(filtered_record.keys()), ':'+', :'.join(filtered_record.keys()))
        c.execute(sql_query, filtered_record)
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')
        raise Exception("Tried to insert in table '{}' record:\n{}\n".format(table, record))

    # if just inserting, return last row ID, which is non-zero if insertion succeeded
    if mode == 'insert':
        return c.lastrowid
    # try to return existing record ID
    elif mode == 'insert-or-ignore':
        if not c.rowcount: # cannot test here on c.lastrowid because of a bug; sometimes when record is ignored a large non-zero value is returned as of sqlite3 version 2.6.0 for python 2.7
            records = get_records(conn, table, 'id', filtered_record)
            if len(records) > 1:
                raise Exception("Cannot retrieve record rowid because query returned multiple records.")
            elif not len(records):
                raise Exception("Tried to insert in table '{}' record:\n{}\n".format(table, record))
            else:
                return records[0]['id']
        else:
            return c.lastrowid
    else:
        raise ValueError("Invalid record addition mode '{}'.\n".format(mode))

def get_records(conn, table, select_fields, record):
    """
    Obtains the record ID.

    Parameters
    ----------
    conn :
        DB connection.

    table : str
        DB table name.

    select_fields : str or list of str
        Name of fields to return.

    record : dict
        Record used to filter DB entries

    """
    if type(select_fields) is str:
        select_fields = [select_fields]

    try:
        c = conn.cursor()
        # build query
        sql_query = "SELECT "+','.join(select_fields)+" FROM %s WHERE "%table
        filters = []
        for rec_field in record:
            filters.append(rec_field+"=:"+rec_field)
        sql_query += ' AND '.join(filters)

        c.execute(sql_query, record)
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')

    return c.fetchall()

def update_records(conn, table, select_filter, update):
    """
    Updates existing records.

    Parameters
    ----------
    conn :
        DB connection.

    table : str
        DB table name.

    select_filter : dict
        Record selection filter.

    update : dict
        Fields to update.
    """
    if not update:
        return
    try:
        c = conn.cursor()
        query_filters = []
        for f in select_filter:
            if type(select_filter[f]) in [str, unicode]:
                sf = '"{}"'.format(select_filter[f])
            elif select_filter[f] is None:
                sf = 'NULL'
            else:
                sf = str(select_filter[f])
            query_filters.append(f+"="+sf)
        update_fields = []
        for f in update:
            if type(update[f]) in [str, unicode]:
                upd = '"{}"'.format(update[f])
            elif update[f] is None:
                upd = 'NULL'
            else:
                upd = str(update[f])
            update_fields.append(f+"="+upd)

        if query_filters:
            sql_query = "UPDATE {} SET {} WHERE {}".format(table, ','.join(update_fields),','.join(query_filters))
        else:
            sql_query = "UPDATE {} SET {}".format(table, ','.join(update_fields))
        c.execute(sql_query)
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')

def sql_query_records(conn, query, output_type = 'dict', parse_dates = None):
    """
    Querry DB records using SQL syntax.
    
    Parameters
    ----------
    conn :
        DB connection.
        
    query : str
        SQL query.

    output_type : str
        Coose from 'dict' or 'dataframe'.

    parse_dates : None or list of str
        If output_type == 'dataframe', convert strings to dates for given column headers.

    Returns
    -------
    1. If output_type == 'dict':
        list of dict
            DB records with query selected fields.

    2. If output_type == 'dataframe':
        pandas.DataFrame
    """
    
    if output_type == 'dict':
        try:
            c = conn.cursor()
            c.execute(query)
            out = [dict(row) for row in c.fetchall()]
        except sqlite3.Error as e:
            util.clrd_print("DB error: {}\n".format(e), 'error')
            out = []
    elif output_type == 'dataframe':
        try:
            out = pd.read_sql_query(query, conn, parse_dates = None)
        except sqlite3.Error as e:
            util.clrd_print("DB error: {}\n".format(e), 'error')
            out = pd.DataFrame()

    return out

def load_aligned_sig(scan_path):
    """
    Loads aligned signals from provided folder. The relative path of the aligned signals file must be:
    "./analysis/aligned_signals.pkl"

    Parameters
    ----------
    scan_path : str
        Folder path as stored in the DB.

    Returns
    -------
    dict or None
        If aligned signals file was loaded successfully, returns dict, None otherwise.
    """
    fpath = os.path.join(scan_path, 'analysis/aligned_signals.pkl')
    sig = None
    if os.path.exists(fpath):
        with open(fpath, 'rb') as f:    
            sig = pkl.load(f)
    else:
        util.clrd_print("Signals file {} is missing.\n".format(fpath), 'error')

    return sig

def save_aligned_sig(sig, scan_path):
    """
    Saves aligned signals file in the selected folder. Note the aligned signals file will be placed in:
    scan_path/analysis/aligned_signals.pkl

    Parameters
    ----------
    sig : dict
        Aligned signals dict.

    scan_path : str
        Folder path as stored in the DB.

    Returns
    -------
    None
    """
    analysis_folder_path = os.path.join(scan_path, 'analysis')
    if not os.path.exists(analysis_folder_path):
        os.mkdir(analysis_folder_path)
    aligned_sig_fpath = os.path.join(analysis_folder_path, 'aligned_signals.pkl')    
    # save to .pkl
    with open(aligned_sig_fpath, 'wb') as f:
        pkl.dump(sig, f)

def load_scan_paths(conn, db_filters = {}):
    """
    Loads all valid (non-excluded) scan paths organized for each cell.

    Parameters
    ----------
    conn : sqlite3.Connection
        DB connection.

    db_filters : dict or str
        DB dataset filters.
        If dict, optional keys are:        

            session_dates : list of str
                Filter for session dates as yyyy-mm-dd

            mouseIDs : list of str
                Filter for mice.

            cellIDs : list of str
                Filter for cell IDs.

            scanIDs : list of int
                Scan numbers.
        If str:
            SQL conditionals that will be appended to the WHERE block.

    Returns
    -------
    dict of list of str
        Scan paths, with dict keys as cell IDs and values as list of scan path strings.
    """
    sql_query = """  SELECT cells.cellID AS cells_cellID,
                            scans.fpath AS scans_fpath
                     FROM   scans, ROIfields, sessions, cells
                     WHERE  scans.ROIfield_id = ROIfields.id AND
                            ROIfields.session_id = sessions.id AND
                            sessions.EP_id = cells.EP_id AND
                            NOT scans.exclude
                """

    # parse DB search filters
    if isinstance(db_filters, dict):
        if 'cellIDs' in db_filters and db_filters['cellIDs']:
            sql_query += " AND cells.cellID IN (%s)"%(','.join(["'"+cid+"'" for cid in db_filters['cellIDs']]))
        if 'mouseIDs' in db_filters and db_filters['mouseIDs']:
            sql_query += " AND mice.mouse_name IN (%s)"%(','.join(["'"+mid+"'" for mid in db_filters['mouseIDs']]))
        if 'scanIDs' in db_filters and db_filters['scanIDs']:
            sql_query += " AND scans.scan IN (%s)"%','.join([str(sid) for sid in db_filters['scanIDs']])
        if 'session_dates' in db_filters and db_filters['session_dates']:
            sql_query += " AND date(sessions.session_datetime) IN (%s)"%(','.join(['date("%s")'%d for d in db_filters['session_dates'].split(',')]))
    elif isinstance(db_filters, basestring):
        sql_query += " "+db_filters
    else:
        raise ValueError

    scan_paths = {}
    # get scan records
    records = sql_query_records(conn, sql_query)
    if not records:
        util.clrd_print("Error: no scan records exists or could be retrieved using filters.\n", 'error')
    for r in records:
        if r['cells_cellID'] not in scan_paths:
            scan_paths[r['cells_cellID']] = []
        scan_paths[r['cells_cellID']].append(r['scans_fpath'])

    return scan_paths

def load_target_signals(conn, targets):
    """
    Loads aligned signals data from ROItargets.
    
    Parameters
    ----------
    conn : sqlite3.Connection
        DB connection.
    
    targets : list of int
        List of ROItarget DB id.

    Returns
    -------
    dict
        Aligned signals with ROItarget.id keys and aligned signal values.
    """
    signals = {}
    for t in targets:
        try:
            c = conn.cursor()
            c.execute("SELECT fpath FROM scans INNER JOIN ROItargets on scans.id = ROItargets.scan_id WHERE ROItargets.id = ?", (t,))        
            record = c.fetchone()
            if record:
                scan_path = record[0]
                fpath = os.path.join(scan_path, 'analysis/aligned_signals.pkl')
                if os.path.exists(fpath):
                    with open(fpath, 'rb') as f:    
                        signals[t] = pkl.load(f)
                else:
                    util.clrd_print("Signals file {} is missing.\n".format(fpath), 'error')
            else:
                util.clrd_print("Cannot load signals. Query for ROItarget ID {} did not return any record.\n".format(t), 'error')
        except sqlite3.Error as e:
            util.clrd_print("DB error: {}\n".format(e), 'error')
    return signals

def get_scans(conn, session_dates = [], mouseIDs = [], cellIDs = [], scanIDs = [], output_type = 'dataframe'):
    """
    Obtains a list of scans with additional information.
    Parameters
    ----------
    conn : sqlite3.Connection
        DB connection.
        
    session_dates : list of str
        Filter for session dates as yyyy-mm-dd

    mouseIDs : list of str
        Filter for mice.

    cellIDs : list of str
        Filter for cell IDs.

    scanIDs : list of int
        Scan numbers.
    
    output_type : str
        Output type: Choose from 'dict' and 'dataframe'.

    Returns
    -------
    list of dict or pandas.DataFrame
    """
    sql_scan_records_query = """  SELECT cells.cellID AS cells_cellID,
                                         sessions.session_date AS sessions_session_date,  
                                         scans.scan AS scans_scan,
                                         scans.fpath AS scans_fpath,
                                         scans.exclude AS scans_exclude
                                  FROM   scans, ROIfields, sessions, cells
                                  WHERE  scans.ROIfield_id = ROIfields.id AND
                                         ROIfields.session_id = sessions.id AND
                                         sessions.EP_id = cells.EP_id{};
                             """

    # parse DB search filters
    search_filters = []
    if cellIDs:
        search_filters.append("cells.cellID IN (%s)"%(','.join(["'"+cid+"'" for cid in cellIDs.split(',')])))
    if session_dates:
        search_filters.append('date(sessions.session_datetime) IN (%s)'%(','.join(['date("%s")'%d for d in session_dates.split(',')])))
    if scanIDs is not None:
        search_filters.append("scans.scan IN (%s)"%scanIDs)
    if mouseIDs:
        search_filters.append("mice.mouse_name IN (%s)"%(','.join(["'"+mid+"'" for mid in mouseIDs])))
    
    if len(search_filters):
        search_filters.insert(0, "")

    # apply search filters to DB query    
    filt_sql_scan_records_query = sql_scan_records_query.format(" AND ".join(search_filters))

    return sql_query_records(conn, filt_sql_scan_records_query, output_type = output_type, parse_dates = ['sessions_session_date'])

def get_ROI_targets(conn, db_filters = {}, selector = [], output_type = 'dataframe'):
    """
    Obtains all ROI scan targets that belong to the same location on the cell, sorted in order of their acquisition time
    and retrieves additional information such as:
        - morphology data if available

    Parameters
    ----------
    conn : sqlite3.Connection
        DB connection.
        
    db_filters : dict or str
        DB dataset filters.
        If dict, optional keys are:        

            session_dates : list of str
                Filter for session dates as yyyy-mm-dd

            mouseIDs : list of str
                Filter for mice.

            cellIDs : list of str
                Filter for cell IDs.

            scanIDs : list of int
                Scan numbers.
        If str:
            SQL conditionals that will be appended to the WHERE block.

    selector : list of str
        Domain specific search of records that contain additional information related to the ROI group. Records that cannot
        be linked are ignored. The following domain selectors are implemented:
        'cell_morph': Complements ROI target information with dendrite morphology information. Adds following fields:
            'cell_morph_section' : str
            'cell_morph_segment' : str
            'cell_morph_segment_class_lvl1' : str
            'cell_morph_segment_class_lvl2' : str
            'cell_morph_segment_dist_to_sec_start' : float
            'cell_morph_section_length' : float
            'cell_morph_soma_distance' : float
            'cell_morph_soma_trunk_location' : float
            'cell_morph_dist_to_soma_trunk_axis' : float
            'cell_morph_soma_trunk_stem_segment': str

    output_type : str
        Output type: Choose from 'dict' and 'dataframe'.

    Returns
    -------
    1. If output_type == 'dict':
        dict with <mouse_name>/<cellID> keys and values:
            dict with with <scanfield ROI>/<ROI target group>/<ROI target name> keys and list of dict with keys and values:
                'ROItargets_id' : int
                    DB ROItargets.id
                'scans_scan' : int
                    Scan #.
                'scans_fpath' : str
                    Folder containing aligned signals file.
                'scans_scan_start_tstamp' : datetime.datetime
                    Scan start timestamp.
                'target_types_type' : str
                    Scanned target type, see valid types in target_types DB table.

    2. If output_type == 'dataframe':
        pandas.DataFrame with column headers returned by SQL query
    """
    # add here more general fields if needed
    select_fields = ["imaging_rigs.rigID",
                     "mice.mouse_name",
                     "cells.cellID",
                     "cells.morph_fpath",
                     "ROIfields.scanfield_ROI",
                     "ROItargets.target_group_name",
                     "ROItargets.target_name",
                     "ROItargets.id",
                     "scans.scan",
                     "scans.scan_start_tstamp",
                     "scans.fpath",
                     "scans.id",
                     "scans.fs",
                     "target_types.type"]

    extra_tables = ['']
    # WARNING: when adding extra tables always add more contraints otherwise many more records will be added
    for s in selector:
        # add here domain specific field selectors and extra tables needed
        if s == 'cell_morph':
            select_fields.extend(["cell_morph.section",
                                  "cell_morph.segment",
                                  "cell_morph.segment_loc",
                                  "cell_morph.segment_class_lvl1",
                                  "cell_morph.segment_class_lvl2",
                                  "cell_morph.segment_dist_to_sec_start",
                                  "cell_morph.section_length",
                                  "cell_morph.soma_distance",
                                  "cell_morph.soma_trunk_location",
                                  "cell_morph.dist_to_soma_trunk_axis",
                                  "cell_morph.soma_trunk_stem_segment"])
            extra_tables.append('cell_morph')
        else:
            raise ValueError
    
    sql_query = "SELECT "+','.join([f+' AS '+f.replace('.','_') for f in select_fields])+\
                """ FROM   mice, EP, sessions, imaging_rigs, ROIfields, scans, ROItargets, cells, target_types {}   
                    WHERE  mice.id = EP.mouse_id AND
                           EP.id = sessions.EP_id AND
                           sessions.id = ROIfields.session_id AND
                           sessions.imaging_rig_id = imaging_rigs.id AND
                           ROIfields.id = scans.ROIfield_id AND
                           scans.id = ROItargets.scan_id AND
                           ROItargets.cell_id = cells.id AND
                           ROItargets.target_type_id = target_types.id AND
                           NOT scans.exclude AND
                           target_types.type NOT IN ("background", "undefined")
                """.format(','.join(extra_tables))

    # add here domain specific conditions to link tables
    for s in selector:
        if s == 'cell_morph':
            sql_query += " AND cell_morph.ROItarget_id = ROItargets.id"
        else:
            raise ValueError

    if isinstance(db_filters, dict):
        if 'cellIDs' in db_filters and db_filters['cellIDs']:
            sql_query += " AND cells.cellID IN (%s)"%(','.join(["'"+cid+"'" for cid in db_filters['cellIDs']]))
        if 'mouseIDs' in db_filters and db_filters['mouseIDs']:
            sql_query += " AND mice.mouse_name IN (%s)"%(','.join(["'"+mid+"'" for mid in db_filters['mouseIDs']]))
        if 'scanIDs' in db_filters and db_filters['scanIDs']:
            sql_query += " AND scans.scan IN (%s)"%','.join([str(sid) for sid in db_filters['scanIDs']])
        if 'session_dates' in db_filters and db_filters['session_dates']:
            sql_query += " AND date(sessions.session_datetime) IN (%s)"%(','.join(['date("%s")'%d for d in db_filters['session_dates'].split(',')]))
    elif isinstance(db_filters, basestring):
        sql_query += " "+db_filters
    else:
        raise ValueError
                  
    # slower than dataframe
    if output_type == 'dict':
        records = []
        try:
            c = conn.cursor()

            #import timeit
            #s = timeit.default_timer()
            c.execute(sql_query)
            query_result = c.fetchall()
            #print(timeit.default_timer()-s)
            #raise Exception
            records = [dict(row) for row in query_result]

        except sqlite3.Error as e:
            util.clrd_print("DB error: {}\n".format(e), 'error')
            return {}

        out = {}
        for r in records:
            mouse_cell = '/'.join([r['mice_mouse_name'], r['cells_cellID']])
            if mouse_cell not in out:
                out[mouse_cell] = {}
            scan_target = '/'.join([r['ROIfields_scanfield_ROI'], r['ROItargets_target_group_name'], r['ROItargets_target_name']])
            if scan_target not in out[mouse_cell]:
                out[mouse_cell][scan_target] = []
            # modify select record fields to make processing easier later
            r['scans_scan_start_tstamp'] = dateutil.parser.parse(r['scans_scan_start_tstamp'])
            # add record
            out[mouse_cell][scan_target].append({k:v for k,v in iteritems(r) if k not in ['mice_mouse_name', 'cells_cellID', 'ROIfields_scanfield_ROI', 
                     'ROItargets_target_group_name', 'ROItargets_target_name']})

        # sort scan groups according to their recording timestamps in ascending order   
        """
        for _,v1 in iteritems(out):
            for _,v2 in iteritems(v1):
                v2.sort(key = lambda x: x['scans_scan_start_tstamp'])
        """
    elif output_type == 'dataframe':
        try:

            #import timeit
            #s = timeit.default_timer()
            out = pd.read_sql_query(sql_query, conn, parse_dates = ['scans_scan_start_tstamp'])
            #print(timeit.default_timer()-s)
            #raise Exception

        except sqlite3.Error as e:
            util.clrd_print("DB error: {}\n".format(e), 'error')
            return pd.DataFrame()
    else:
        raise ValueError()

    return out

def get_ROItarget_morph_info(db_conn, ROItarget_id):
    """
    Obtains morphology info for a given ROItarget id in the DB.

    Parameters
    ----------
    db_conn : sqlite3.Connection
        DB connection.

    ROItarget_id : int
        ROI target ID from ROItargets table.

    Returns
    -------
    pandas.DataFraame of 1 row
    """
    sql_query = """ SELECT 
                            cells.cellID AS cells_cellID,
                            ROIfields.scanfield_ROI AS ROIfields_scanfield_ROI,
                            scans.scan AS scans_scan,
                            
                            target_types.type AS target_types_type,
                            ROItargets.target_group_name AS ROItargets_target_group_name,
                            ROItargets.target_name AS ROItargets_target_name,

                            cell_morph.section AS cell_morph_section,
                            cell_morph.segment_dist_to_sec_start AS cell_morph_segment_dist_to_sec_start,
                            cell_morph.section_length AS cell_morph_section_length,
                            cell_morph.soma_distance AS cell_morph_soma_distance,
                            cell_morph.soma_trunk_location AS cell_morph_soma_trunk_location,
                            cell_morph.dist_to_soma_trunk_axis AS cell_morph_dist_to_soma_trunk_axis,
                            cell_morph.soma_trunk_stem_segment AS cell_morph_soma_trunk_stem_segment,
                            
                            mice.mouse_name AS mice_mouse_name,
                            EP.EP_datetime AS EP_EP_datetime,
                            sessions.session_date AS sessions_session_date
                            
                    FROM    ROItargets, ROIfields, cell_morph, target_types, scans, sessions, EP, mice, cells
                    WHERE   ROItargets.id = %d AND
                            cell_morph.ROItarget_id = ROItargets.id AND
                            ROItargets.target_type_id = target_types.id AND
                            ROItargets.scan_id = scans.id AND
                            scans.ROIfield_id = ROIfields.id AND 
                            ROIfields.session_id = sessions.id AND
                            sessions.EP_id = EP.id AND 
                            EP.mouse_id = mice.id AND
                            cells.EP_id = EP.id AND
                            ROItargets.cell_id = cells.id
                """%ROItarget_id

    try:
        query_df = pd.read_sql_query(sql_query, db_conn)
        if len(query_df) > 1: # must return no more than one row from DB
            util.clrd_print("DB error: ROItarget ID %d query returned more than one record."%ROItarget_id, 'error')
            print(query_df)
            raise Exception()
        return query_df
    except sqlite3.Error as e:
        util.clrd_print("DB error: {}\n".format(e), 'error')
        return pd.Series()
        
                
    