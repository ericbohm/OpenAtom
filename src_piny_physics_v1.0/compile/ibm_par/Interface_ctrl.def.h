/* DEFS: readonly MDINTEGRATE readonly_mdintegrate;
 */
extern MDINTEGRATE readonly_mdintegrate;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_readonly_mdintegrate(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|readonly_mdintegrate;
}
#endif

/* DEFS: readonly MDATOMS readonly_mdatoms;
 */
extern MDATOMS readonly_mdatoms;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_readonly_mdatoms(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|readonly_mdatoms;
}
#endif

/* DEFS: readonly MDINTER readonly_mdinter;
 */
extern MDINTER readonly_mdinter;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_readonly_mdinter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|readonly_mdinter;
}
#endif

/* DEFS: readonly MDINTRA readonly_mdintra;
 */
extern MDINTRA readonly_mdintra;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_readonly_mdintra(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|readonly_mdintra;
}
#endif

/* DEFS: readonly GENERAL_DATA readonly_general_data;
 */
extern GENERAL_DATA readonly_general_data;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_readonly_general_data(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|readonly_general_data;
}
#endif

/* DEFS: readonly CP readonly_cp;
 */
extern CP readonly_cp;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_readonly_cp(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|readonly_cp;
}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registerInterface_ctrl(void)
{
  static int _done = 0; if(_done) return; _done = 1;
  CkRegisterReadonly("readonly_mdintegrate","MDINTEGRATE",sizeof(readonly_mdintegrate),(void *) &readonly_mdintegrate,__xlater_roPup_readonly_mdintegrate);

  CkRegisterReadonly("readonly_mdatoms","MDATOMS",sizeof(readonly_mdatoms),(void *) &readonly_mdatoms,__xlater_roPup_readonly_mdatoms);

  CkRegisterReadonly("readonly_mdinter","MDINTER",sizeof(readonly_mdinter),(void *) &readonly_mdinter,__xlater_roPup_readonly_mdinter);

  CkRegisterReadonly("readonly_mdintra","MDINTRA",sizeof(readonly_mdintra),(void *) &readonly_mdintra,__xlater_roPup_readonly_mdintra);

  CkRegisterReadonly("readonly_general_data","GENERAL_DATA",sizeof(readonly_general_data),(void *) &readonly_general_data,__xlater_roPup_readonly_general_data);

  CkRegisterReadonly("readonly_cp","CP",sizeof(readonly_cp),(void *) &readonly_cp,__xlater_roPup_readonly_cp);

}
#endif
