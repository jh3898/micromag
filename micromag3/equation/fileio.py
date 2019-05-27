
import os
import numpy as np
import re

class DataSaver:
    comment_symbol = '# '

    def __init__(self, sim, filename, entities = None):
        self.sim = sim
        self.filename = filename

        precision = 12
        charwidth = 18
        self.float_format = '%' + str(charwidth) + '.' + str(precision) + 'g '
        self.string_format = '%' + str(charwidth) + 's '


        self.entities = {
                    'step':{'unit': '<>',
                            'get': lambda sim: sim.driver.step,
                            'header': 'step'},
                    'time' :{
                             'unit': '<>',
                             'get': lambda sim : sim.driver.t,
                             'header': 'time'
                             },
                    'm': {'unit': '<>',
                            'get': lambda sim: sim.compute_average(),
                            'header': ('m_x', 'm_y', 'm_z')}
                }
        if entities is not None:
            self.entities = entities

        self.save_head = False
        self.entity_order = self.default_entity_order()

    def default_entity_order(self):
        keys = sorted(self.entities.keys())

        if 'time' in keys:
            keys.remove('time')
            return ['time'] + sorted(keys)

        elif 'step' in keys:
            keys.remove('step')
            return ['step'] + sorted(keys)
        return list(keys)
    def update_entity_order(self):
        self.entity_order = self.default_entity_order()

    def headers(self):
        """
        return line one and two of ndt data file as string
        :return: 
        """###
        line1 = [self.comment_symbol]
        line2 = [self.comment_symbol]
        for entityname in self.entity_order:
            headers = self.entities[entityname]['header']

            if isinstance(headers, str):
                headers = [headers]
            for head in headers:
                line1.append(self.string_format % head)
                line2.append(self.string_format % self.entities[entityname]['unit'])

        return "".join(line1) + '\n' + ''.join(line2) + '\n'


    def save(self):
        if not self.save_head:
            f = open(self.filename, 'w')
            f.write(self.headers())
            f.close()
            self.save_head = True

        with open(self.filename, 'a') as f:
            f.write(' ' * len(self.comment_symbol))

            for entityname in self.entity_order:
                value = self.entities[entityname]['get'](self.sim)
                if isinstance(value, np.ndarray):
                    for v in value:
                        f.write(self.float_format % v)
                elif isinstance(value, float) or isinstance(value, int):
                    f.write(self.float_format % value)
                elif value is None:
                    f.write(self.string_format % 'nan')
                else:
                    msg = 'Can only deal with numpy arrays, float and int' + \
                        ', but type is %s' % type(value)
                    raise NotImplementedError(msg)


            f.write('\n')


if __name__ == '__main__':
    pass







