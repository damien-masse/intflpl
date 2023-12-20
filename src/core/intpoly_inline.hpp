/**
 * Polyhedra using intervals (include file)
 * ----------------------------------------------------------------------------
 *  \date       2023
 *  \author     Damien Massé
 *  \copyright  Copyright 2023
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */


inline int IntPoly::get_dim() const { return this->dim; }
inline int IntPoly::size() const { return this->dim; }
inline int IntPoly::get_nbcsts() const { return this->csts.size(); }
inline int IntPoly::get_not_flat_dim() const { return this->dim_not_flat; }
inline const IntervalVector &IntPoly::getBox() const { return this->Box; }
inline const IntervalVector &IntPoly::box() const { return this->Box; }
inline const CstrVectMap &IntPoly::getCsts() const { return this->csts; }
inline bool IntPoly::is_empty() const { return pstatus[INTPOLY_EMPTY]; }
inline bool IntPoly::is_flat() const { return this->dim_not_flat<this->dim; }
inline bool IntPoly::is_bounded() const {return !pstatus[INTPOLY_UNBOUNDED]; }
inline bool IntPoly::is_unbounded() const {return pstatus[INTPOLY_UNBOUNDED]; }
inline bool IntPoly::is_box() const { return this->csts.size()==0; }
inline const Interval& IntPoly::operator[](unsigned int i) const { return this->Box[i]; }
inline Vector IntPoly::mid() const { return this->Box.mid(); }
inline void IntPoly::set_empty() {
        if (this->pstatus[INTPOLY_INVALID]) return;
        this->pstatus=1<<INTPOLY_EMPTY;
        this->dim_not_flat=-1;
        this->Inc.set_empty();
        this->Box.set_empty(); this->csts.clear(); this->minimized=0; }
