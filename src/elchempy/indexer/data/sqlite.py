"""
This module contains the sql interface for data manipulation.
"""

import array
import functools
import json
import logging
logger = logging.getLogger(__name__)

import sqlite3

import pandas as pd

from elchempy.indexer.data import DATABASE




def






def with_connection(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        # If cursor exists we just move on
        if kwargs.get('cursor'):
            return func(*args, **kwargs)

        db_path = kwargs.get('db_path', DATABASE)
        db_path = db_path if db_path else DATABASE
        conn = sqlite3.connect(
            db_path
        )  # TODO remove 'str' call when dropping P3.6
        conn.row_factory = sqlite3.Row

        try:
            # Get a cursor object
            cursor = conn.cursor()
            cursor.execute('PRAGMA foreign_keys = ON')
            ret = func(*args, **kwargs, cursor=cursor)

        except sqlite3.IntegrityError as e:
            conn.rollback()
            raise ParsingError(e)

        except sqlite3.InterfaceError as e:
            conn.rollback()
            raise ParsingError(e)

        else:
            conn.commit()

        finally:
            conn.close()

        return ret

    return wrapper


# ---------------------- General functions


def _upload_one_all_columns(
    cursor,
    table_name,
    table_id,
    columns,
    input_dict,
    overwrite,
    print_string,
    verbose,
    **kwargs,
):
    """Insert or overwrite a list of things in a table."""
    to_insert = [table_id] + columns

    if overwrite:
        sql_com = build_update(
            table=table_name, to_set=columns, where=[table_id]
        )
    else:
        sql_com = build_insert(table=table_name, to_insert=to_insert)

    # Upload or modify data
    insert_dict = {key: input_dict.get(key) for key in to_insert}
    try:
        cursor.execute(sql_com, insert_dict)
    except sqlite3.Error as e:
        raise type(e)(
            f"Error inserting dict {insert_dict}. Original error:\n {e}"
        )

    if verbose:
        # Print success
        logger.info(f"{print_string} uploaded {insert_dict.get(table_id)}")


def _get_all_no_id(
    cursor,
    table_name,
    table_id,
    print_string,
    verbose,
    **kwargs,
):
    """Get all elements from a table as a dictionary, excluding id."""
    try:
        cursor.execute("""SELECT * FROM """ + table_name)
    except sqlite3.Error as e:
        raise type(e)(
            f"Error getting data from {table_name}. Original error:\n {e}"
        )

    values = []
    for row in cursor:
        val = dict(zip(row.keys(), row))
        val.pop(table_id)
        values.append(val)

    return values


def _delete_by_id(
    cursor,
    table_name,
    table_id,
    element_id,
    print_string,
    verbose,
    **kwargs,
):
    """Delete elements in a table by using their ID."""

    # Check if exists
    ids = cursor.execute(
        build_select(table=table_name, to_select=[table_id], where=[table_id]),
        {
            table_id: element_id
        }
    ).fetchone()

    if ids is None:
        raise sqlite3.IntegrityError(
            f"Element to delete ({element_id}) does not exist in table {table_name}."
        )

    try:
        cursor.execute(
            build_delete(table=table_name, where=[table_id]),
            {table_id: element_id}
        )
    except sqlite3.Error as e:
        raise type(e)(
            f"Error deleting {element_id} from {table_name}. Original error:\n {e}"
        )

    if verbose:
        # Print success
        logger.info(f"Success, deleted {print_string}: {element_id}")


# ---------------------- Adsorbates


@with_connection
def index_to_db(
    adsorbate, db_path=None, overwrite=False, verbose=True, **kwargs
):
    """
    Upload an adsorbate to the database.
    If overwrite is set to true, the adsorbate is overwritten.
    Overwrite is done based on adsorbate.name
    Parameters
    ----------
    adsorbate : Adsorbate
        Adsorbate class to upload to the database.
    db_path : str, None
        Path to the database. If none is specified, internal database is used.
    overwrite : bool
        Whether to upload the adsorbate or overwrite it.
        WARNING: Overwrite is done on ALL fields.
    verbose : bool
        Extra information printed to console.
    """
    cursor = kwargs['cursor']

    # If we need to overwrite, we find the id of existing adsorbate.
    if overwrite:
        ids = cursor.execute(
            build_select(table='index', to_select=['id'], where=['name']),
            {
                'name': adsorbate.name
            }
        ).fetchone()
        if ids is None:
            raise sqlite3.IntegrityError(
                f"Adsorbate to overwrite ({adsorbate.name}) does not exist in database."
            )
        ads_id = ids[0]
    # If overwrite is not specified, we upload it to the adsorbates table
    else:
        cursor.execute(
            build_insert(table="adsorbates", to_insert=['name']),
            {'name': adsorbate.name}
        )
        ads_id = cursor.lastrowid

    # pload or modify data in the associated tablesU
    properties = adsorbate.to_dict()
    del properties['name']  # no need for this

    if properties:
        if overwrite:
            # Delete existing properties
            _delete_by_id(
                cursor,
                'adsorbate_properties',
                'ads_id',
                ads_id,
                'adsorbate properties',
                verbose,
            )

        for prop, val in properties.items():

            sql_insert = build_insert(
                table='adsorbate_properties',
                to_insert=['ads_id', 'type', 'value'],
            )

            if not isinstance(val, (list, set, tuple)):
                val = [val]

            for vl in val:
                try:
                    cursor.execute(
                        sql_insert, {
                            'ads_id': ads_id,
                            'type': prop,
                            'value': vl
                        }
                    )
                except sqlite3.InterfaceError as e:
                    raise type(e)(
                        f"Cannot process property {prop}: {vl}"
                        f"Original error:\n{e}"
                    )

    # Add to existing list
    if overwrite:
        if adsorbate in ADSORBATE_LIST:
            ADSORBATE_LIST.remove(adsorbate.name)
    ADSORBATE_LIST.append(adsorbate)

    if verbose:
        # Print success
        logger.info(f"Adsorbate uploaded: '{adsorbate.name}'")
