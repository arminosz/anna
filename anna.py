import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, Crippen
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image, ImageTk
import csv
import os


class AminoAcidAnalogGenerator:

    def __init__(self, root):
        self.root = root
        self.root.title('ANNA')
        self.root.geometry('1400x800')
        self.amino_acids = {'Alanina': 'N[C@@H](C)C(=O)O', 'Valina':
            'N[C@@H](C(C)C)C(=O)O', 'Leucina': 'N[C@@H](CC(C)C)C(=O)O',
            'Isoleucina': 'N[C@@H]([C@@H](C)CC)C(=O)O', 'Serina':
            'N[C@@H](CO)C(=O)O', 'Treonina': 'N[C@@H]([C@H](O)C)C(=O)O',
            'Cisteína': 'N[C@@H](CS)C(=O)O', 'Metionina':
            'N[C@@H](CCSC)C(=O)O', 'Fenilalanina':
            'N[C@@H](Cc1ccccc1)C(=O)O', 'Tirosina':
            'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'Triptofano':
            'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O', 'Histidina':
            'N[C@@H](Cc1c[nH]cn1)C(=O)O', 'Asparagina':
            'N[C@@H](CC(=O)N)C(=O)O', 'Glutamina':
            'N[C@@H](CCC(=O)N)C(=O)O', 'Ácido Aspártico':
            'N[C@@H](CC(=O)O)C(=O)O', 'Ácido Glutâmico':
            'N[C@@H](CCC(=O)O)C(=O)O', 'Lisina': 'N[C@@H](CCCCN)C(=O)O',
            'Arginina': 'N[C@@H](CCCNC(=N)N)C(=O)O', 'Prolina':
            'N1[C@@H](CCC1)C(=O)O', 'Glicina': 'NCC(=O)O'}
        self.modifications = {'Alifática': {'Metilação': ['[CH3]',
            '[CH2][CH3]', '[CH2][CH2][CH3]'], 'Ramificação': [
            '[CH]([CH3])[CH3]', '[C]([CH3])([CH3])[CH3]'], 'Ciclopropil': [
            '[CH2]C1CC1']}, 'Aromática': {'Substituição benzílica': [
            'c1ccc(F)cc1', 'c1ccc(Cl)cc1', 'c1ccc(Br)cc1'],
            'Heteroaromática': ['c1ccncc1', 'c1cccnc1', 'c1ccoc1'],
            'Naftaleno': ['c1ccc2ccccc2c1']}, 'Hidroxilação': {
            'Hidroxil primário': ['[CH2]O', '[CH2][CH2]O'],
            'Hidroxil secundário': ['[CH](O)[CH3]', '[CH](O)[CH2][CH3]'],
            'Diol': ['[CH](O)[CH2]O', '[CH](O)[CH](O)[CH3]']}, 'Metilação':
            {'N-metilação': ['[NH][CH3]', '[NH]([CH3])[CH3]'],
            'O-metilação': ['O[CH3]'], 'Metil terminal': ['[CH2][CH3]',
            '[CH2][CH2][CH3]']}, 'Halogenação': {'Fluoração': ['F',
            '[CH2]F', '[CHF2]'], 'Cloração': ['Cl', '[CH2]Cl'], 'Bromação':
            ['Br', '[CH2]Br']}, 'Grupos funcionais': {'Nitro': [
            '[N+](=O)[O-]'], 'Ciano': ['C#N'], 'Carbonil': ['C(=O)[CH3]',
            'C(=O)O'], 'Amida': ['C(=O)N', 'C(=O)N([CH3])[CH3]'], 'Éster':
            ['C(=O)O[CH3]', 'C(=O)O[CH2][CH3]']}}
        self.current_analogs = []
        self.setup_gui()

    def setup_gui(self):
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.generation_frame = ttk.Frame(notebook)
        notebook.add(self.generation_frame, text='Geração de Análogos')
        self.visualization_frame = ttk.Frame(notebook)
        notebook.add(self.visualization_frame, text=
            'Visualização e Propriedades')
        self.setup_generation_tab()
        self.setup_visualization_tab()

    def setup_generation_tab(self):
        left_frame = ttk.Frame(self.generation_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        ttk.Label(left_frame, text='Aminoácido Base:', font=('Arial', 12,
            'bold')).pack(anchor=tk.W, pady=(0, 5))
        self.aa_var = tk.StringVar(value='Alanina')
        aa_combo = ttk.Combobox(left_frame, textvariable=self.aa_var,
            values=list(self.amino_acids.keys()), state='readonly', width=25)
        aa_combo.pack(anchor=tk.W, pady=(0, 15))
        ttk.Label(left_frame, text='Tipo de Modificação:', font=('Arial', 
            12, 'bold')).pack(anchor=tk.W, pady=(0, 5))
        self.mod_var = tk.StringVar(value='Alifática')
        mod_combo = ttk.Combobox(left_frame, textvariable=self.mod_var,
            values=list(self.modifications.keys()), state='readonly', width=25)
        mod_combo.pack(anchor=tk.W, pady=(0, 15))
        mod_combo.bind('<<ComboboxSelected>>', self.update_submodifications)
        ttk.Label(left_frame, text='Submodificação:', font=('Arial', 12,
            'bold')).pack(anchor=tk.W, pady=(0, 5))
        self.submod_var = tk.StringVar()
        self.submod_combo = ttk.Combobox(left_frame, textvariable=self.
            submod_var, state='readonly', width=25)
        self.submod_combo.pack(anchor=tk.W, pady=(0, 15))
        self.update_submodifications()
        ttk.Label(left_frame, text='Máximo de Análogos:', font=('Arial', 12,
            'bold')).pack(anchor=tk.W, pady=(0, 5))
        self.max_analogs_var = tk.IntVar(value=10)
        max_spin = tk.Spinbox(left_frame, from_=1, to=20, textvariable=self
            .max_analogs_var, width=23)
        max_spin.pack(anchor=tk.W, pady=(0, 15))
        generate_btn = ttk.Button(left_frame, text='Gerar Análogos',
            command=self.generate_analogs, style='Accent.TButton')
        generate_btn.pack(pady=20, fill=tk.X)
        ttk.Label(left_frame, text='Filtros:', font=('Arial', 12, 'bold')
            ).pack(anchor=tk.W, pady=(20, 5))
        filter_frame = ttk.LabelFrame(left_frame, text='Filtro por logP')
        filter_frame.pack(fill=tk.X, pady=5)
        ttk.Label(filter_frame, text='Min:').grid(row=0, column=0, padx=5,
            pady=5)
        self.logp_min_var = tk.DoubleVar(value=-5.0)
        tk.Spinbox(filter_frame, from_=-10, to=10, increment=0.1,
            textvariable=self.logp_min_var, width=8).grid(row=0, column=1,
            padx=5, pady=5)
        ttk.Label(filter_frame, text='Max:').grid(row=0, column=2, padx=5,
            pady=5)
        self.logp_max_var = tk.DoubleVar(value=5.0)
        tk.Spinbox(filter_frame, from_=-10, to=10, increment=0.1,
            textvariable=self.logp_max_var, width=8).grid(row=0, column=3,
            padx=5, pady=5)
        filter_frame2 = ttk.LabelFrame(left_frame, text='Filtro por Massa (Da)'
            )
        filter_frame2.pack(fill=tk.X, pady=5)
        ttk.Label(filter_frame2, text='Min:').grid(row=0, column=0, padx=5,
            pady=5)
        self.mw_min_var = tk.DoubleVar(value=50.0)
        tk.Spinbox(filter_frame2, from_=50, to=1000, increment=10,
            textvariable=self.mw_min_var, width=8).grid(row=0, column=1,
            padx=5, pady=5)
        ttk.Label(filter_frame2, text='Max:').grid(row=0, column=2, padx=5,
            pady=5)
        self.mw_max_var = tk.DoubleVar(value=500.0)
        tk.Spinbox(filter_frame2, from_=50, to=1000, increment=10,
            textvariable=self.mw_max_var, width=8).grid(row=0, column=3,
            padx=5, pady=5)
        ttk.Button(left_frame, text='Aplicar Filtros', command=self.
            apply_filters).pack(pady=10, fill=tk.X)
        export_frame = ttk.LabelFrame(left_frame, text='Exportação')
        export_frame.pack(fill=tk.X, pady=10)
        ttk.Button(export_frame, text='Exportar para CSV', command=self.
            export_csv).pack(pady=5, fill=tk.X)
        ttk.Button(export_frame, text='Exportar para SDF', command=self.
            export_sdf).pack(pady=5, fill=tk.X)
        right_frame = ttk.Frame(self.generation_frame)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10,
            pady=10)
        ttk.Label(right_frame, text='Análogos Gerados:', font=('Arial', 14,
            'bold')).pack(anchor=tk.W)
        columns = 'ID', 'SMILES', 'logP', 'MW', 'HBD', 'HBA', 'PSA'
        self.tree = ttk.Treeview(right_frame, columns=columns, show=
            'headings', height=20)
        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=100 if col != 'SMILES' else 200)
        scrollbar = ttk.Scrollbar(right_frame, orient=tk.VERTICAL, command=
            self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.bind('<<TreeviewSelect>>', self.on_analog_select)

    def setup_visualization_tab(self):
        structure_frame = ttk.LabelFrame(self.visualization_frame, text=
            'Estrutura 2D')
        structure_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=
            10, pady=10)
        self.structure_label = ttk.Label(structure_frame, text=
            'Selecione um análogo para visualizar')
        self.structure_label.pack(expand=True)
        self.copy_frame = ttk.Frame(structure_frame)
        self.copy_frame.pack(pady=10)
        self.smiles_var = tk.StringVar()
        smiles_entry = ttk.Entry(self.copy_frame, textvariable=self.
            smiles_var, width=40, state='readonly')
        smiles_entry.pack(side=tk.LEFT, padx=5)
        ttk.Button(self.copy_frame, text='Copiar SMILES', command=self.
            copy_smiles).pack(side=tk.LEFT)
        props_frame = ttk.LabelFrame(self.visualization_frame, text=
            'Propriedades Fisicoquímicas')
        props_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=10, pady=10)
        self.props_tree = ttk.Treeview(props_frame, columns=('Propriedade',
            'Valor'), show='headings', height=15)
        self.props_tree.heading('Propriedade', text='Propriedade')
        self.props_tree.heading('Valor', text='Valor')
        self.props_tree.column('Propriedade', width=200)
        self.props_tree.column('Valor', width=150)
        self.props_tree.pack(fill=tk.BOTH, expand=True)

    def update_submodifications(self, event=None):
        mod_type = self.mod_var.get()
        if mod_type in self.modifications:
            submods = list(self.modifications[mod_type].keys())
            self.submod_combo['values'] = submods
            if submods:
                self.submod_var.set(submods[0])

    def generate_analogs(self):
        try:
            aa_name = self.aa_var.get()
            mod_type = self.mod_var.get()
            submod = self.submod_var.get()
            max_analogs = self.max_analogs_var.get()
            if not all([aa_name, mod_type, submod]):
                messagebox.showerror('Erro',
                    'Por favor, selecione todas as opções.')
                return
            base_smiles = self.amino_acids[aa_name]
            base_mol = Chem.MolFromSmiles(base_smiles)
            if base_mol is None:
                messagebox.showerror('Erro',
                    'SMILES do aminoácido base inválido.')
                return
            analogs = []
            modifications_list = self.modifications[mod_type][submod]
            for i, modification in enumerate(modifications_list[:max_analogs]):
                analog_smiles = self.apply_modification(base_smiles,
                    modification, aa_name)
                if analog_smiles:
                    mol = Chem.MolFromSmiles(analog_smiles)
                    if mol:
                        properties = self.calculate_properties(mol)
                        analog_data = {'ID': f'{aa_name}_{submod}_{i + 1}',
                            'SMILES': analog_smiles, 'Molecule': mol, **
                            properties}
                        analogs.append(analog_data)
            self.current_analogs = analogs
            self.update_analog_table()
            messagebox.showinfo('Sucesso',
                f'{len(analogs)} análogos gerados com sucesso!')
        except Exception as e:
            messagebox.showerror('Erro',
                f'Erro na geração de análogos: {str(e)}')

    def apply_modification(self, base_smiles, modification, aa_name):
        try:
            base_mol = Chem.MolFromSmiles(base_smiles)
            if aa_name in ['Fenilalanina', 'Tirosina', 'Triptofano']:
                if 'c1ccc' in modification:
                    if aa_name == 'Fenilalanina':
                        return f'N[C@@H](C{modification})C(=O)O'
                    elif aa_name == 'Tirosina':
                        return (
                            f"N[C@@H](Cc1ccc(O)c({modification.replace('c1ccc', '')})c1)C(=O)O"
                            )
            elif aa_name in ['Alanina', 'Valina', 'Leucina']:
                if '[CH3]' in modification:
                    return base_smiles.replace('C)C(=O)O',
                        f'C{modification})C(=O)O')
            if aa_name == 'Alanina':
                return f'N[C@@H](C{modification})C(=O)O'
            elif aa_name == 'Serina':
                return f'N[C@@H](C{modification}O)C(=O)O'
            elif aa_name == 'Fenilalanina':
                return f'N[C@@H](Cc1ccc({modification})cc1)C(=O)O'
            elif aa_name == 'Lisina':
                return f'N[C@@H](CCCC{modification}N)C(=O)O'
            else:
                return base_smiles.replace(')C(=O)O', f'{modification})C(=O)O')
        except Exception as e:
            print(f'Erro na modificação: {e}')
            return None

    def calculate_properties(self, mol):
        try:
            properties = {'logP': round(Crippen.MolLogP(mol), 2), 'MW':
                round(Descriptors.MolWt(mol), 2), 'HBD': Descriptors.
                NumHDonors(mol), 'HBA': Descriptors.NumHAcceptors(mol),
                'PSA': round(Descriptors.TPSA(mol), 2), 'RotBonds':
                Descriptors.NumRotatableBonds(mol), 'Rings': Descriptors.
                RingCount(mol), 'Charge': self.estimate_charge_at_ph(mol, 
                7.4), 'Solubility': self.estimate_solubility(mol)}
            return properties
        except Exception as e:
            print(f'Erro no cálculo de propriedades: {e}')
            return {}

    def estimate_charge_at_ph(self, mol, ph):
        basic_patterns = [Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]'), Chem.
            MolFromSmarts('[NX3;H0;!$(NC=O)]')]
        acidic_patterns = [Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]'),
            Chem.MolFromSmarts('[#16X2H1]')]
        basic_count = sum(len(mol.GetSubstructMatches(pattern)) for pattern in
            basic_patterns if pattern)
        acidic_count = sum(len(mol.GetSubstructMatches(pattern)) for
            pattern in acidic_patterns if pattern)
        charge = basic_count - acidic_count
        return round(charge, 1)

    def estimate_solubility(self, mol):
        logp = Crippen.MolLogP(mol)
        psa = Descriptors.TPSA(mol)
        if logp < 0 and psa > 60:
            return 'Alta'
        elif logp < 2 and psa > 40:
            return 'Moderada'
        elif logp < 4:
            return 'Baixa'
        else:
            return 'Muito Baixa'

    def update_analog_table(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        for analog in self.current_analogs:
            self.tree.insert('', tk.END, values=(analog['ID'], analog[
                'SMILES'][:50] + '...' if len(analog['SMILES']) > 50 else
                analog['SMILES'], analog.get('logP', 'N/A'), analog.get(
                'MW', 'N/A'), analog.get('HBD', 'N/A'), analog.get('HBA',
                'N/A'), analog.get('PSA', 'N/A')))

    def apply_filters(self):
        if not self.current_analogs:
            return
        logp_min = self.logp_min_var.get()
        logp_max = self.logp_max_var.get()
        mw_min = self.mw_min_var.get()
        mw_max = self.mw_max_var.get()
        filtered_analogs = []
        for analog in self.current_analogs:
            logp = analog.get('logP', 0)
            mw = analog.get('MW', 0)
            if logp_min <= logp <= logp_max and mw_min <= mw <= mw_max:
                filtered_analogs.append(analog)
        original_analogs = self.current_analogs.copy()
        self.current_analogs = filtered_analogs
        self.update_analog_table()
        messagebox.showinfo('Filtros',
            f'{len(filtered_analogs)} análogos passaram pelos filtros.')

    def on_analog_select(self, event):
        selection = self.tree.selection()
        if not selection:
            return
        item = self.tree.item(selection[0])
        analog_id = item['values'][0]
        selected_analog = None
        for analog in self.current_analogs:
            if analog['ID'] == analog_id:
                selected_analog = analog
                break
        if selected_analog:
            self.display_analog_details(selected_analog)

    def display_analog_details(self, analog):
        self.smiles_var.set(analog['SMILES'])
        try:
            mol = analog['Molecule']
            img = Draw.MolToImage(mol, size=(400, 300))
            img_tk = ImageTk.PhotoImage(img)
            self.structure_label.configure(image=img_tk, text='')
            self.structure_label.image = img_tk
        except Exception as e:
            self.structure_label.configure(text=
                f'Erro na visualização: {str(e)}', image='')
        for item in self.props_tree.get_children():
            self.props_tree.delete(item)
        properties = [('ID', analog['ID']), ('SMILES', analog['SMILES']), (
            'logP', analog.get('logP', 'N/A')), ('Massa Molecular (Da)',
            analog.get('MW', 'N/A')), ('Doadores de H', analog.get('HBD',
            'N/A')), ('Aceitadores de H', analog.get('HBA', 'N/A')), (
            'PSA (Ų)', analog.get('PSA', 'N/A')), ('Ligações Rotáveis',
            analog.get('RotBonds', 'N/A')), ('Número de Anéis', analog.get(
            'Rings', 'N/A')), ('Carga estimada (pH 7.4)', analog.get(
            'Charge', 'N/A')), ('Solubilidade estimada', analog.get(
            'Solubility', 'N/A'))]
        for prop, value in properties:
            self.props_tree.insert('', tk.END, values=(prop, value))

    def copy_smiles(self):
        smiles = self.smiles_var.get()
        if smiles:
            self.root.clipboard_clear()
            self.root.clipboard_append(smiles)
            messagebox.showinfo('Copiado',
                'SMILES copiado para área de transferência!')

    def export_csv(self):
        if not self.current_analogs:
            messagebox.showwarning('Aviso', 'Nenhum análogo para exportar.')
            return
        filename = filedialog.asksaveasfilename(defaultextension='.csv',
            filetypes=[('CSV files', '*.csv'), ('All files', '*.*')])
        if filename:
            try:
                data = []
                for analog in self.current_analogs:
                    row = {'ID': analog['ID'], 'SMILES': analog['SMILES'],
                        'logP': analog.get('logP', ''), 'MW': analog.get(
                        'MW', ''), 'HBD': analog.get('HBD', ''), 'HBA':
                        analog.get('HBA', ''), 'PSA': analog.get('PSA', ''),
                        'RotBonds': analog.get('RotBonds', ''), 'Rings':
                        analog.get('Rings', ''), 'Charge_pH7.4': analog.get
                        ('Charge', ''), 'Solubility': analog.get(
                        'Solubility', '')}
                    data.append(row)
                df = pd.DataFrame(data)
                df.to_csv(filename, index=False)
                messagebox.showinfo('Sucesso',
                    f'Dados exportados para {filename}')
            except Exception as e:
                messagebox.showerror('Erro', f'Erro na exportação: {str(e)}')

    def export_sdf(self):
        if not self.current_analogs:
            messagebox.showwarning('Aviso', 'Nenhum análogo para exportar.')
            return
        filename = filedialog.asksaveasfilename(defaultextension='.sdf',
            filetypes=[('SDF files', '*.sdf'), ('All files', '*.*')])
        if filename:
            try:
                writer = Chem.SDWriter(filename)
                for analog in self.current_analogs:
                    mol = analog['Molecule']
                    for prop, value in analog.items():
                        if prop not in ['Molecule', 'SMILES']:
                            mol.SetProp(str(prop), str(value))
                    writer.write(mol)
                writer.close()
                messagebox.showinfo('Sucesso',
                    f'Estruturas exportadas para {filename}')
            except Exception as e:
                messagebox.showerror('Erro',
                    f'Erro na exportação SDF: {str(e)}')


def main():
    root = tk.Tk()
    style = ttk.Style()
    style.theme_use('clam')
    style.configure('Accent.TButton', background='#0078d4', foreground='white')
    style.map('Accent.TButton', background=[('active', '#106ebe')])
    app = AminoAcidAnalogGenerator(root)
    root.update_idletasks()
    x = root.winfo_screenwidth() // 2 - root.winfo_width() // 2
    y = root.winfo_screenheight() // 2 - root.winfo_height() // 2
    root.geometry(f'+{x}+{y}')
    root.mainloop()


if __name__ == '__main__':
    main()