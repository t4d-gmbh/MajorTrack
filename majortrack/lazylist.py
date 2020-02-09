from __future__ import absolute_import, print_function, unicode_literals
import os
from datetime import datetime
from collections import MutableSequence
import pickle
import shutil
import itertools
import atexit


class lazy_list(MutableSequence):
    def __init__(
            self,
            lazy_iterator=itertools.count(start=0, step=1),
            lazy_name='{0}'.format(
                    datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
                    ),
            temp_path=os.path.join(os.getcwd(), '.lazy_list_objs'),
            *args):
        self._lazy_iterator = lazy_iterator
        self._name = lazy_name
        self._temp_path = os.path.join(temp_path, self._name)
        if os.path.exists(self._temp_path):
            free = itertools.count(start=2, step=1)
            while os.path.exists(self._temp_path):
                _name = '{0}_{1}'.format(self._name, next(free))
                self._temp_path = os.path.join(temp_path, _name)
            self._name = _name
        os.makedirs(self._temp_path)
        atexit.register(self._cleanup)
        self._list = list()
        self.extend(list(args))

    @property
    def _lazy_next(self):
        return os.path.join(
                self._temp_path,
                '{0}.p'.format(next(self._lazy_iterator))
                )

    def __len__(self):
        return len(self._list)

    def __getitem__(self, ii):
        with open(self._list[ii], 'rb') as fobj:
            return pickle.load(fobj, encoding='bytes')

    def __delitem__(self, ii):
        os.remove(self._list[ii])
        del self._list[ii]

    def __setitem__(self, ii, value):
        _file_name = self._lazy_next
        assert os.path.isfile(_file_name) is False
        with open(_file_name, 'wb') as fobj:
            pickle.dump(value, fobj)
        del value
        self._list[ii] = _file_name

    def insert(self, ii, value):
        _file_name = self._lazy_next
        assert os.path.isfile(_file_name) is False
        with open(_file_name, 'wb') as fobj:
            pickle.dump(value, fobj)
        del value
        self._list.insert(ii, _file_name)

    def __str__(self):
        return str(self._list)

    def _cleanup(self,):
        shutil.rmtree(self._temp_path)
